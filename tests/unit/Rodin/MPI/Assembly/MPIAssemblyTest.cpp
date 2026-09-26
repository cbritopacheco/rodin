/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for the MPI assembly module.
 *
 * These tests verify that the MPI assembly infrastructure works correctly:
 * - MPIIteration iterates over local mesh elements
 * - Assembly::MPI<IndexMap<...>, DirichletBC<...>> assembles boundary conditions
 *   on a distributed mesh
 * - For 1-rank runs the assembled boundary DOF set is consistent with the
 *   sequential assembly
 * - P1, H1 orders 1–4, P0g, and real/complex/vector traces reach DOF owners
 *
 * Run with mpirun -n 1/2/3/4/8 as registered in CMakeLists.txt.
 */
#include <algorithm>
#include <numeric>
#include <set>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/Assembly.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/MPI/Variational/P0/P0.h>
#include <Rodin/MPI/Variational/H1/H1.h>
#include <Rodin/MPI/Variational/P0g/P0g.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------------------------------------------------------
// Global MPI handles (initialized in main())
// ---------------------------------------------------------------------------
static boost::mpi::environment* g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /**
   * @brief Creates a local mesh with all incidences required for sharding
   * and boundary-DOF assembly.
   */
  static Mesh<Context::Local> makeShardableMesh(
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    if (D > 2)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
    }
    return mesh;
  }

  /**
   * @brief Distributes a uniform-grid mesh from rank 0.
   */
  static Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx,
    Polytope::Type type, std::initializer_list<size_t> shape,
    std::vector<Index>* physicalBoundary = nullptr)
  {
    const auto& comm = ctx.getCommunicator();

    Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(type, shape);
      if (physicalBoundary)
        for (auto face = localMesh.getBoundary(); face; ++face)
          physicalBoundary->push_back(face->getIndex());
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);
    if (physicalBoundary)
      boost::mpi::broadcast(comm, *physicalBoundary, 0);
    return mesh;
  }

  template <class FES>
  static std::set<Index> requiredDOFs(const FES& fes)
  {
    std::set<Index> result;
    Index begin, end;
    fes.getOwnershipRange(begin, end);
    for (Index i = begin; i < end; ++i)
      result.insert(i);
    const auto& mesh = fes.getMesh();
    for (auto cell = mesh.getCell(); cell; ++cell)
      if (mesh.getShard().isOwned(mesh.getDimension(), cell->getIndex()))
        for (Index dof : fes.getDOFs(mesh.getDimension(), cell->getIndex()))
          result.insert(dof);
    return result;
  }

  static const char* polytopeName(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Tetrahedron:
        return "Tetrahedron";
      case Polytope::Type::Hexahedron:
        return "Hexahedron";
      case Polytope::Type::Pyramid:
        return "Pyramid";
      case Polytope::Type::Wedge:
        return "Wedge";
      case Polytope::Type::Segment:
        return "Segment";
      case Polytope::Type::Triangle:
        return "Triangle";
      case Polytope::Type::Quadrilateral:
        return "Quadrilateral";
      default:
        return "Other";
    }
  }
}

namespace Rodin::Tests::Unit
{
  class MPITraceGeometryTest : public testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::MPI> makeMesh(
        const Context::MPI& ctx, std::vector<Index>* physicalBoundary = nullptr) const
      {
        const auto geometry = GetParam();
        const size_t dim = Polytope::Traits(geometry).getDimension();
        if (dim == 1)
          return distributeFromRoot(ctx, geometry, {17}, physicalBoundary);
        if (dim == 2)
          return distributeFromRoot(ctx, geometry, {5, 5}, physicalBoundary);
        const size_t n = geometry == Polytope::Type::Tetrahedron ? 9 : 5;
        return distributeFromRoot(ctx, geometry, {n, n, n}, physicalBoundary);
      }
  };

  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPITraceGeometryTest,
    testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) { return polytopeName(info.param); });

  /**
   * Every supported trace DOF must be constrained on its owning rank.
   * On this three-rank tetrahedral partition, face-local assembly alone
   * omitted respectively 2, 2, 4, 6, and 8 owned DOFs for P1 and H1 orders
   * one through four. Vector spaces duplicate the omissions per component.
   */
  TEST_P(MPITraceGeometryTest, BoundaryConstraintsReachOwnersAcrossSpaces)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = makeMesh(ctx);
    auto probe = [&](const auto& fes, const auto& value, const char* name) {
      TrialFunction u(fes);
      DirichletBC dbc(u, value);
      dbc.assemble();
      using Scalar = typename std::remove_cvref_t<decltype(fes)>::ScalarType;
      const auto& corrected = std::get<IndexMap<Scalar>>(dbc.getDOFs());
      std::set<Index> localBoundary;
      for (auto it = mesh.getBoundary(); it; ++it)
      {
        const auto& indices = fes.getDOFs(mesh.getDimension() - 1, it->getIndex());
        localBoundary.insert(indices.begin(), indices.end());
      }
      std::vector<Index> localVector(localBoundary.begin(), localBoundary.end());
      std::vector<std::vector<Index>> gathered;
      boost::mpi::all_gather(world, localVector, gathered);
      std::set<Index> boundary;
      for (const auto& indices : gathered)
        boundary.insert(indices.begin(), indices.end());
      Index begin = 0, end = 0;
      fes.getOwnershipRange(begin, end);
      int missing = 0;
      for (const Index index : boundary)
      {
        if (begin <= index && index < end)
          missing += !corrected.contains(index);
      }
      const int total = boost::mpi::all_reduce(world, missing, std::plus<int>());
      EXPECT_FALSE(boundary.empty()) << name;
      EXPECT_EQ(total, 0) << name;
    };
    const RealFunction realValue(1);
    const ComplexFunction complexValue(Complex(1, 2));
    const VectorFunction vectorValue{realValue, realValue, realValue};
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    probe(p1, realValue, "P1");
    H1<1, Real, Mesh<Context::MPI>> h1(std::integral_constant<size_t, 1>{}, mesh);
    probe(h1, realValue, "H1<1>");
    H1<2, Real, Mesh<Context::MPI>> h2(std::integral_constant<size_t, 2>{}, mesh);
    probe(h2, realValue, "H1<2>");
    H1<3, Real, Mesh<Context::MPI>> h3(std::integral_constant<size_t, 3>{}, mesh);
    probe(h3, realValue, "H1<3>");
    H1<4, Real, Mesh<Context::MPI>> h4(std::integral_constant<size_t, 4>{}, mesh);
    probe(h4, realValue, "H1<4>");
    P0g<Real, Mesh<Context::MPI>> p0g(mesh);
    probe(p0g, realValue, "P0g");
    P0g<Complex, Mesh<Context::MPI>> complexP0g(mesh);
    probe(complexP0g, complexValue, "P0g complex");
    H1<2, Complex, Mesh<Context::MPI>> complexH2(
      std::integral_constant<size_t, 2>{}, mesh);
    probe(complexH2, complexValue, "H1<2> complex");
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP1(mesh, 3);
    probe(vectorP1, vectorValue, "P1 vector");
    H1<2, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH2(
      std::integral_constant<size_t, 2>{}, mesh, 3);
    probe(vectorH2, vectorValue, "H1<2> vector");
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP0g(mesh, 3);
    probe(vectorP0g, vectorValue, "P0g vector");
    P1<Complex, Mesh<Context::MPI>> complexP1(mesh);
    probe(complexP1, complexValue, "P1 complex");
    H1<1, Complex, Mesh<Context::MPI>> complexH1(
      std::integral_constant<size_t, 1>{}, mesh);
    probe(complexH1, complexValue, "H1<1> complex");
    H1<1, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH1(
      std::integral_constant<size_t, 1>{}, mesh, 3);
    probe(vectorH1, vectorValue, "H1<1> vector");
    H1<3, Complex, Mesh<Context::MPI>> complexH3(
      std::integral_constant<size_t, 3>{}, mesh);
    probe(complexH3, complexValue, "H1<3> complex");
    H1<3, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH3(
      std::integral_constant<size_t, 3>{}, mesh, 3);
    probe(vectorH3, vectorValue, "H1<3> vector");
    H1<4, Complex, Mesh<Context::MPI>> complexH4(
      std::integral_constant<size_t, 4>{}, mesh);
    probe(complexH4, complexValue, "H1<4> complex");
    H1<4, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH4(
      std::integral_constant<size_t, 4>{}, mesh, 3);
    probe(vectorH4, vectorValue, "H1<4> vector");
  }

  /**
   * Every rank owning a slave or assembling a cell using it needs its row.
   * A three-rank tetrahedral P2 partition has boundary slaves on ranks
   * other than their DOF owner; owner-only reconciliation is insufficient
   * because an off-process slave can occur in a locally assembled cell.
   * The same invariant is checked for real, complex, and vector traces.
   */
  TEST_P(MPITraceGeometryTest, LocalTraceAssemblyIsNoncollective)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = makeMesh(ctx);
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    H1<2, Real, Mesh<Context::MPI>> p2(std::integral_constant<size_t, 2>{}, mesh);
    auto probe = [&](const auto& fes) {
      // TrialFunction allocates a distributed GridFunction collectively.
      // Only constraint assembly, not distributed-field construction, is local.
      TrialFunction u(fes);
      auto value = DirichletBC(u, RealFunction(2));
      auto affine = DirichletBC(u, RealFunction(-1) * u, RealFunction(4));
      if (g_world->rank() == 0)
      {
        value.assemble();
        affine.assemble();
        const auto& prescribed =
          std::get<DirichletBCBase<Real>::ValueDOFs>(value.getDOFs());
        EXPECT_FALSE(prescribed.empty());
        EXPECT_EQ(prescribed.size(), affine.getIdentificationValues().size());
      }
      g_world->barrier();
    };
    probe(p1);
    probe(p2);
  }

  /** Owned-cell extraction must not constrain partition interfaces (formerly 10 extra P1 DOFs). */
  TEST_P(MPITraceGeometryTest, SubMeshBoundaryConstraintsMatchPhysicalBoundary)
  {
    Context::MPI ctx(*g_env, *g_world);
    std::vector<Index> boundary;
    auto parent = GetParam() == Polytope::Type::Triangle
      ? distributeFromRoot(ctx, GetParam(), {9, 9}, &boundary)
      : makeMesh(ctx, &boundary);
    const std::set<Index> physical(boundary.begin(), boundary.end());
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(parent);
    const size_t D = parent.getDimension();
    for (auto cell = parent.getCell(); cell; ++cell)
      if (parent.getShard().isOwned(D, cell->getIndex()))
        builder.include(D, cell->getIndex());
    auto sub = builder.finalize();
    const Mesh<Context::MPI>& mesh = sub;
    const auto probe = [&](const auto& fes) {
      TrialFunction u(fes);
      auto dbc = DirichletBC(u, RealFunction(1));
      dbc.assemble();
      const auto& values = std::get<IndexMap<Real>>(dbc.getDOFs());
      const auto required = requiredDOFs(fes);
      std::set<Index> expected;
      for (auto face = mesh.getFace(); face; ++face)
        if (physical.contains(
              mesh.getShard().getPolytopeMap(D - 1).left.at(face->getIndex())))
          for (Index dof : fes.getDOFs(D - 1, face->getIndex()))
            if (required.contains(dof))
              expected.insert(dof);
      EXPECT_EQ(values.size(), expected.size());
      for (Index dof : expected)
        EXPECT_TRUE(values.contains(dof));
    };
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    H1<2, Real, Mesh<Context::MPI>> p2(std::integral_constant<size_t, 2>{}, mesh);
    probe(p1);
    probe(p2);
  }

  TEST_P(MPITraceGeometryTest, IdentificationRowsReachRequiredDOFs)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = makeMesh(ctx);
    auto probe = [&](const auto& fes, const char* name) {
      TrialFunction u(fes);
      TrialFunction v(fes);
      using Space = std::remove_cvref_t<decltype(fes)>;
      const auto defect = [] {
        if constexpr (std::is_same_v<typename Space::RangeType,
                        Math::SpatialVector<Real>>)
          return VectorFunction{RealFunction(2), RealFunction(3), RealFunction(4)};
        else if constexpr (std::is_same_v<typename Space::ScalarType, Complex>)
          return ComplexFunction(Complex(2, 3));
        else
          return RealFunction(2);
      }();
      DirichletBC dbc(u, -v, defect);
      dbc.assemble();
      using Scalar = typename std::remove_cvref_t<decltype(fes)>::ScalarType;
      using IdentifiedDOFs = typename DirichletBCBase<Scalar>::IdentifiedDOFs;
      const auto& rows = std::get<IdentifiedDOFs>(dbc.getDOFs());
      std::vector<Index> local;
      // Independent oracle: enumerate the original owned boundary faces.
      for (auto it = mesh.getBoundary(); it; ++it)
        for (Index dof : fes.getDOFs(mesh.getDimension() - 1, it->getIndex()))
          local.push_back(dof);
      std::vector<std::vector<Index>> gathered;
      boost::mpi::all_gather(world, local, gathered);
      std::set<Index> all;
      for (const auto& keys : gathered)
        all.insert(keys.begin(), keys.end());
      ASSERT_FALSE(all.empty()) << name;
      const auto required = requiredDOFs(fes);
      const auto& offsets = dbc.getIdentificationValues();
      EXPECT_EQ(offsets.size(), rows.size());
      for (const Index slave : all)
      {
        const auto it = rows.find(slave);
        if (!required.contains(slave))
        {
          EXPECT_EQ(it, rows.end()) << name << " unrelated slave=" << slave;
          continue;
        }
        ASSERT_NE(it, rows.end())
          << name << " rank=" << world.rank() << " slave=" << slave;
        EXPECT_TRUE(offsets.contains(slave)) << name;
        EXPECT_GT(it->second.first.size(), 0) << name;
        EXPECT_EQ(it->second.first.size(), it->second.second.size()) << name;
      }
    };
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    probe(p1, "P1 real");
    H1<1, Real, Mesh<Context::MPI>> h1(std::integral_constant<size_t, 1>{}, mesh);
    probe(h1, "H1<1> real");
    H1<2, Real, Mesh<Context::MPI>> h2(std::integral_constant<size_t, 2>{}, mesh);
    probe(h2, "H1<2> real");
    H1<3, Real, Mesh<Context::MPI>> h3(std::integral_constant<size_t, 3>{}, mesh);
    probe(h3, "H1<3> real");
    H1<4, Real, Mesh<Context::MPI>> h4(std::integral_constant<size_t, 4>{}, mesh);
    probe(h4, "H1<4> real");
    P0g<Real, Mesh<Context::MPI>> p0g(mesh);
    probe(p0g, "P0g real");
    P0g<Complex, Mesh<Context::MPI>> complexP0g(mesh);
    probe(complexP0g, "P0g complex");
    H1<2, Complex, Mesh<Context::MPI>> complexH2(
      std::integral_constant<size_t, 2>{}, mesh);
    probe(complexH2, "H1<2> complex");
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP1(mesh, 3);
    probe(vectorP1, "P1 vector");
    H1<2, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH2(
      std::integral_constant<size_t, 2>{}, mesh, 3);
    probe(vectorH2, "H1<2> vector");
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP0g(mesh, 3);
    probe(vectorP0g, "P0g vector");
    P1<Complex, Mesh<Context::MPI>> complexP1(mesh);
    probe(complexP1, "P1 complex");
    H1<1, Complex, Mesh<Context::MPI>> complexH1(
      std::integral_constant<size_t, 1>{}, mesh);
    probe(complexH1, "H1<1> complex");
    H1<1, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH1(
      std::integral_constant<size_t, 1>{}, mesh, 3);
    probe(vectorH1, "H1<1> vector");
    H1<3, Complex, Mesh<Context::MPI>> complexH3(
      std::integral_constant<size_t, 3>{}, mesh);
    probe(complexH3, "H1<3> complex");
    H1<3, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH3(
      std::integral_constant<size_t, 3>{}, mesh, 3);
    probe(vectorH3, "H1<3> vector");
    H1<4, Complex, Mesh<Context::MPI>> complexH4(
      std::integral_constant<size_t, 4>{}, mesh);
    probe(complexH4, "H1<4> complex");
    H1<4, Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorH4(
      std::integral_constant<size_t, 4>{}, mesh, 3);
    probe(vectorH4, "H1<4> vector");
  }

  /** Certifies the mesh metadata independently of constraint assembly. */
  TEST_P(MPITraceGeometryTest, HaloContainsRequiredPhysicalBoundaryIncidence)
  {
    Context::MPI ctx(*g_env, *g_world);
    std::vector<Index> parentBoundary;
    auto mesh = makeMesh(ctx, &parentBoundary);
    const std::set<Index> physical(parentBoundary.begin(), parentBoundary.end());
    const auto& shard = mesh.getShard();
    const size_t dim = mesh.getDimension();
    const auto probe = [&](const auto& space, const char* name) {
      SCOPED_TRACE(name);
      // No DBC assembler or new exchange participates in this oracle.
      std::vector<std::pair<Index, Index>> incidences;
      IndexMap<std::set<Index>> visible;
      std::set<Index> ownedFaceDOFs, required;
      Index begin, end;
      space.getOwnershipRange(begin, end);
      for (Index dof = begin; dof < end; ++dof)
        required.insert(dof);
      for (auto cell = mesh.getCell(); cell; ++cell)
        if (shard.isOwned(dim, cell->getIndex()))
          for (Index dof : space.getDOFs(dim, cell->getIndex()))
            required.insert(dof);

      size_t falseFaces = 0, relevantFalseFaces = 0;
      for (auto face = mesh.getFace(); face; ++face)
      {
        const Index local = face->getIndex();
        const Index global = shard.getPolytopeMap(dim - 1).left.at(local);
        if (physical.contains(global))
        {
          EXPECT_TRUE(shard.isBoundary(local));
          for (Index dof : space.getDOFs(dim - 1, local))
          {
            visible[dof].insert(global);
            if (shard.isOwned(dim - 1, local))
            {
              incidences.emplace_back(dof, global);
              ownedFaceDOFs.insert(dof);
            }
          }
        }
        else if (shard.isBoundary(local))
        {
          ++falseFaces;
          for (Index dof : space.getDOFs(dim - 1, local))
            relevantFalseFaces += required.contains(dof);
        }
      }
      std::vector<std::vector<std::pair<Index, Index>>> gathered;
      boost::mpi::all_gather(*g_world, incidences, gathered);
      IndexMap<std::set<Index>> expected;
      for (const auto& records : gathered)
        for (const auto& [dof, face] : records)
          expected[dof].insert(face);
      size_t missing = 0, filtered = 0;
      for (Index dof : required)
      {
        const auto it = expected.find(dof);
        if (it == expected.end())
          continue;
        for (Index face : it->second)
          missing += !visible[dof].contains(face);
        filtered += !ownedFaceDOFs.contains(dof);
      }
      EXPECT_EQ(missing, 0u);
      EXPECT_EQ(relevantFalseFaces, 0u);
      const size_t globalMissing =
        boost::mpi::all_reduce(*g_world, missing, std::plus<size_t>());
      const size_t globalFiltered =
        boost::mpi::all_reduce(*g_world, filtered, std::plus<size_t>());
      const size_t globalFalse =
        boost::mpi::all_reduce(*g_world, falseFaces, std::plus<size_t>());
      if (g_world->rank() == 0)
        std::cout << "[halo-audit] " << polytopeName(GetParam()) << ' ' << name
                  << " missing incidences=" << globalMissing
                  << " DOFs missed by owned-face filter=" << globalFiltered
                  << " artificial shard-boundary faces=" << globalFalse << '\n';
    };
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    probe(p1, "P1");
    const auto orders = [&]<size_t K>() {
      H1<K, Real, Mesh<Context::MPI>> scalar(std::integral_constant<size_t, K>{}, mesh);
      probe(scalar, "H1 real");
      H1<K, Complex, Mesh<Context::MPI>> complex(
        std::integral_constant<size_t, K>{}, mesh);
      probe(complex, "H1 complex");
      H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> vector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      probe(vector, "H1 vector");
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
  }

  /** Certifies the mesh metadata independently of constraint assembly. */
  TEST_P(MPITraceGeometryTest, EntityOwnersAndHalosAgree)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = makeMesh(ctx);
    const auto& shard = mesh.getShard();
    const int rank = g_world->rank();
    for (size_t dim = 0; dim <= mesh.getDimension(); ++dim)
    {
      // gid -> (declared owner, ownership flag), independently on each rank.
      using Record = std::pair<Index, std::pair<Index, bool>>;
      std::vector<Record> local;
      for (Index i = 0; i < shard.getPolytopeCount(dim); ++i)
      {
        const bool owned = shard.isOwned(dim, i);
        const Index owner = owned ? static_cast<Index>(rank) : shard.getOwner(dim).at(i);
        local.emplace_back(shard.getPolytopeMap(dim).left.at(i), std::pair{owner, owned});
      }
      std::vector<std::vector<Record>> gathered;
      boost::mpi::all_gather(*g_world, local, gathered);
      IndexMap<Index> owners;
      IndexMap<IndexSet> holders;
      for (size_t peer = 0; peer < gathered.size(); ++peer)
        for (const auto& [gid, declaration] : gathered[peer])
        {
          holders[gid].insert(peer);
          if (declaration.second)
          {
            EXPECT_EQ(declaration.first, peer);
            EXPECT_TRUE(owners.emplace(gid, peer).second)
              << "duplicate owner gid=" << gid;
          }
        }
      for (const auto& records : gathered)
        for (const auto& [gid, declaration] : records)
        {
          const auto owner = owners.find(gid);
          EXPECT_NE(owner, owners.end());
          if (owner != owners.end())
          {
            EXPECT_EQ(owner->second, declaration.first);
          }
        }
      for (Index i = 0; i < shard.getPolytopeCount(dim); ++i)
        if (shard.isOwned(dim, i))
        {
          auto expected = holders.at(shard.getPolytopeMap(dim).left.at(i));
          expected.erase(rank);
          const auto halo = shard.getHalo(dim).find(i);
          if (halo == shard.getHalo(dim).end())
          {
            EXPECT_TRUE(expected.empty());
          }
          else
          {
            EXPECT_EQ(halo->second, expected);
          }
        }
    }
  }

  /** Affine offsets and identification rows must have the same rank scope. */
  TEST(Assembly_MPI_DirichletBC, AffineIdentificationValuesReachRequiredDOFs)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, {9, 9, 9});
    H1<2, Real, Mesh<Context::MPI>> fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);
    DirichletBC dbc(u, -v, RealFunction(2));
    dbc.assemble();
    const auto& rows = std::get<DirichletBCBase<Real>::IdentifiedDOFs>(dbc.getDOFs());
    const auto& values = dbc.getIdentificationValues();
    ASSERT_FALSE(rows.empty());
    for (const auto& [slave, row] : rows)
    {
      const auto it = values.find(slave);
      ASSERT_NE(it, values.end()) << "rank=" << world.rank() << " slave=" << slave;
      EXPECT_NEAR(it->second, 2, 1e-12);
    }

    H1<2, Complex, Mesh<Context::MPI>> complexFES(
      std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction complexU(complexFES);
    TrialFunction complexV(complexFES);
    const Complex prescribed(2, 3);
    DirichletBC complexDBC(complexU, -complexV, ComplexFunction(prescribed));
    complexDBC.assemble();
    const auto& complexRows =
      std::get<DirichletBCBase<Complex>::IdentifiedDOFs>(complexDBC.getDOFs());
    const auto& complexValues = complexDBC.getIdentificationValues();
    ASSERT_FALSE(complexRows.empty());
    for (const auto& [slave, row] : complexRows)
    {
      const auto it = complexValues.find(slave);
      ASSERT_NE(it, complexValues.end())
        << "rank=" << world.rank() << " complex slave=" << slave;
      EXPECT_NEAR(std::abs(it->second - prescribed), 0, 1e-12);
    }
  }

  /** A P0g trace selected only on remote faces still reaches rank 0. */
  TEST(Assembly_MPI_DirichletBC, P0gRemoteBoundaryReachesOwner)
  {
    const auto& world = *g_world;
    if (world.size() < 2)
      GTEST_SKIP() << "Requires a boundary face owned by a nonzero rank.";
    Context::MPI ctx(*g_env, world);
    Sharder<Context::MPI> sharder(ctx);
    constexpr Attribute selected = 97;
    if (world.rank() == 0)
    {
      auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {9, 9});
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(world.size()));
      const size_t dim = localMesh.getDimension();
      sharder.shard(partitioner);
      const auto& rankZeroFaces = sharder.getShards()[0].getPolytopeMap(dim - 1).right;
      for (auto it = localMesh.getBoundary(); it; ++it)
      {
        if (rankZeroFaces.find(it->getIndex()) == rankZeroFaces.end())
        {
          localMesh.setAttribute({dim - 1, it->getIndex()}, selected);
          break;
        }
      }
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mpiMesh = sharder.gather(0);
    size_t localSelected = 0;
    for (auto it = mpiMesh.getBoundary(); it; ++it)
      if (it->getAttribute() == selected)
        ++localSelected;
    size_t visibleSelected = 0;
    for (auto it = mpiMesh.getFace(); it; ++it)
      if (it->getAttribute() == selected)
        ++visibleSelected;
    const size_t selectedCount =
      boost::mpi::all_reduce(world, localSelected, std::plus<size_t>());
    P0g<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);
    auto dbc = DirichletBC(u, RealFunction(2));
    dbc.on(selected);
    dbc.assemble();
    const auto& values = std::get<IndexMap<Real>>(dbc.getDOFs());
    EXPECT_GT(selectedCount, 0u);
    if (world.rank() == 0)
    {
      EXPECT_EQ(localSelected, 0u);
      EXPECT_EQ(visibleSelected, 0u);
      ASSERT_TRUE(values.contains(0));
      EXPECT_EQ(values.at(0), 2);
    }
  }

  /** Tagged interior faces use topological indices; reselection clears old rows. */
  TEST(Assembly_MPI_DirichletBC, TaggedInteriorAndEmptyReselection)
  {
    Context::MPI ctx(*g_env, *g_world);
    Sharder<Context::MPI> sharder(ctx);
    constexpr Attribute selected = 98, absent = 99;
    if (g_world->rank() == 0)
    {
      auto parent = makeShardableMesh(Polytope::Type::Triangle, {9, 9});
      for (auto face = parent.getFace(); face; ++face)
        if (!parent.isBoundary(face->getIndex()))
          parent.setAttribute({1, face->getIndex()}, selected);
      BalancedCompactPartitioner partitioner(parent);
      partitioner.partition(static_cast<size_t>(g_world->size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);
    auto probe = [&](const auto& fes) {
      TrialFunction u(fes);
      auto value = DirichletBC(u, RealFunction(2));
      auto affine = DirichletBC(u, -u, RealFunction(4));
      value.on(selected).assemble();
      affine.on(selected).assemble();
      std::vector<Index> local;
      for (auto face = mesh.getFace(); face; ++face)
        if (mesh.getShard().isOwned(1, face->getIndex()) &&
          face->getAttribute() == selected)
          for (Index dof : fes.getDOFs(1, face->getIndex()))
            local.push_back(dof);
      std::vector<std::vector<Index>> gathered;
      boost::mpi::all_gather(*g_world, local, gathered);
      const auto required = requiredDOFs(fes);
      std::set<Index> expected;
      for (const auto& indices : gathered)
        for (Index dof : indices)
          if (required.contains(dof))
            expected.insert(dof);
      const auto& values = std::get<IndexMap<Real>>(value.getDOFs());
      const auto& rows =
        std::get<DirichletBCBase<Real>::IdentifiedDOFs>(affine.getDOFs());
      EXPECT_EQ(values.size(), expected.size());
      EXPECT_EQ(rows.size(), expected.size());
      EXPECT_EQ(affine.getIdentificationValues().size(), expected.size());
      for (Index dof : expected)
      {
        EXPECT_TRUE(values.contains(dof));
        EXPECT_TRUE(rows.contains(dof));
        EXPECT_TRUE(affine.getIdentificationValues().contains(dof));
      }
      value.on(absent).assemble();
      affine.on(absent).assemble();
      EXPECT_TRUE(values.empty());
      EXPECT_TRUE(rows.empty());
      EXPECT_TRUE(affine.getIdentificationValues().empty());
    };
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    H1<2, Real, Mesh<Context::MPI>> p2(std::integral_constant<size_t, 2>{}, mesh);
    P0g<Real, Mesh<Context::MPI>> p0g(mesh);
    probe(p1);
    probe(p2);
    probe(p0g);
  }

  /** Linear and affine parts use the same source even for varying face data. */
  TEST(Assembly_MPI_DirichletBC, AffinePartsUseTheSameSource)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {9, 9});
    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);
    const RealFunction data([](const Point& p) { return 2 + p(0) + 3 * p(1); });
    auto dbc = DirichletBC(u, data * v, data);
    dbc.assemble();
    const auto& rows = std::get<DirichletBCBase<Real>::IdentifiedDOFs>(dbc.getDOFs());
    const auto& values = dbc.getIdentificationValues();
    ASSERT_EQ(rows.size(), 1u);
    ASSERT_EQ(values.size(), 1u);
    ASSERT_EQ(rows.at(0).first.size(), 1);
    EXPECT_EQ(rows.at(0).first[0], 0u);
    EXPECT_EQ(rows.at(0).second[0], values.at(0));
    std::vector<Real> gathered;
    boost::mpi::all_gather(*g_world, values.at(0), gathered);
    for (Real value : gathered)
      EXPECT_EQ(value, values.at(0));
  }

  /** Small nonzero master coefficients must not be pruned by MPI assembly. */
  TEST(Assembly_MPI_DirichletBC, TinyMasterCoefficientIsRetained)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {9, 9});
    P1<Real, Mesh<Context::MPI>> scalar(mesh);
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vector(mesh, 2);
    TrialFunction u(scalar);
    TrialFunction v(vector);
    const Real tiny = 1e-30;
    auto dbc = DirichletBC(u, v.x() + RealFunction(tiny) * v.y());
    dbc.assemble();
    const auto& rows = std::get<DirichletBCBase<Real>::IdentifiedDOFs>(dbc.getDOFs());
    for (const auto& [slave, row] : rows)
    {
      EXPECT_EQ(row.first.size(), 2);
      bool found = false;
      for (Index j = 0; j < static_cast<Index>(row.second.size()); ++j)
        found |= row.second[j] == tiny;
      EXPECT_TRUE(found) << "slave=" << slave;
    }
  }

  /** Distinct vector dimensions share a getDOFs scratch buffer, not numbering. */
  TEST(Assembly_MPI_DirichletBC, DistinctVectorSpacesPreserveSlaveIndices)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {9, 9});
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> slave(mesh, 2);
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> master(mesh, 3);
    // Collective query: ranks have different numbers of rows below.
    const size_t masterSize = master.getSize();
    TrialFunction u(slave);
    TrialFunction v(master);
    Math::Matrix<Real> projection = Math::Matrix<Real>::Zero(2, 3);
    projection(0, 0) = 1;
    projection(1, 1) = 1;
    auto dbc = DirichletBC(u, MatrixFunction(projection) * v);
    dbc.assemble();
    const auto& rows = std::get<DirichletBCBase<Real>::IdentifiedDOFs>(dbc.getDOFs());
    std::vector<Index> boundary;
    for (auto it = mesh.getBoundary(); it; ++it)
      for (Index dof : slave.getDOFs(1, it->getIndex()))
        boundary.push_back(dof);
    std::vector<std::vector<Index>> gathered;
    boost::mpi::all_gather(*g_world, boundary, gathered);
    const auto local = requiredDOFs(slave);
    for (Index i = 0; i < slave.getShard().getSize(); ++i)
    {
      const Index global = slave.getGlobalIndex(i);
      EXPECT_EQ(slave.getLocalIndex(global), Optional<Index>(i));
    }
    EXPECT_EQ(slave.getLocalIndex(slave.getSize()), std::nullopt);
    std::set<Index> expected;
    for (const auto& indices : gathered)
      for (Index dof : indices)
        if (local.contains(dof))
          expected.insert(dof);
    EXPECT_EQ(rows.size(), expected.size());
    for (Index dof : expected)
      EXPECT_TRUE(rows.contains(dof)) << "slave=" << dof;
    for (const auto& [dof, row] : rows)
    {
      EXPECT_TRUE(expected.contains(dof));
      ASSERT_EQ(row.first.size(), 1);
      EXPECT_LT(row.first[0], masterSize);
      EXPECT_EQ(row.second[0], 1);
    }
  }

  /** An empty master row is still a constraint, not an absent slave. */
  TEST(Assembly_MPI_DirichletBC, ZeroLinearPartRetainsAffineConstraint)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {9, 9});
    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);
    auto dbc = DirichletBC(u, RealFunction(0) * v, RealFunction(2));
    dbc.assemble();
    const auto& rows = std::get<DirichletBCBase<Real>::IdentifiedDOFs>(dbc.getDOFs());
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows.at(0).first.size(), 0);
    EXPECT_EQ(rows.at(0).second.size(), 0);
    EXPECT_EQ(dbc.getIdentificationValues().at(0), 2);
  }
  // =========================================================================
  // MPIIteration
  // =========================================================================

  /**
   * @brief MPIIteration over a Triangle mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, TriangleMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    // At least one cell on every rank for a 4x4 mesh (32 triangles total)
    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief MPIIteration over a Quadrilateral mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, QuadrilateralMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 4, 4 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief Global cell count equals the sum of local counts across all ranks.
   */
  TEST(Assembly_MPI_Iteration, GlobalCellCount_MatchesLocal_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    // Count only owned cells on this rank (ghost cells must be excluded to
    // avoid double-counting when reducing across ranks).
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();
    size_t localCount = 0;
    {
      Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
      for (auto it = iter.getIterator(); it; ++it)
      {
        if (shard.isOwned(D, it->getIndex()))
          ++localCount;
      }
    }

    // Reduce to root
    size_t totalCount = 0;
    boost::mpi::reduce(world, localCount, totalCount, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      // A 4x4 uniform triangle grid has 2*(4-1)*(4-1) = 18 triangles
      EXPECT_EQ(totalCount, 18u);
    }
  }

  // =========================================================================
  // Assembly::MPI<IndexMap<Real>, DirichletBC<...>>
  // =========================================================================

  /**
   * @brief MPI DirichletBC assembly yields at least one boundary DOF entry
   * on a Triangle mesh with a constant boundary condition.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_NonEmpty_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(1.0));
    dbc.assemble();

    // At least globally some DOFs must be fixed on the boundary
    const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();
    size_t globalFixed = 0;
    boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      EXPECT_GT(globalFixed, 0u);
    }
  }

  /**
   * @brief MPI DirichletBC assembly: boundary DOF count is consistent with
   * the sequential assembly on the same mesh for a 1-rank run.
   */
  TEST(Assembly_MPI_DirichletBC, SingleRank_MatchesSequentialCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P1<Real, Mesh<Context::MPI>> mpiFes(mpiMesh);
    TrialFunction uMPI(mpiFes);
    DirichletBC dbcMPI(uMPI, RealFunction(1.0));
    dbcMPI.assemble();
    const size_t mpiFixed = std::get<IndexMap<Real>>(dbcMPI.getDOFs()).size();

    // Sequential reference
    auto localMesh = makeShardableMesh(Polytope::Type::Triangle, { 4, 4 });
    P1 seqFes(localMesh);
    TrialFunction uSeq(seqFes);
    DirichletBC dbcSeq(uSeq, RealFunction(1.0));
    dbcSeq.assemble();
    const size_t seqFixed = std::get<IndexMap<Real>>(dbcSeq.getDOFs()).size();

    EXPECT_EQ(mpiFixed, seqFixed);
  }

  /**
   * @brief MPI DirichletBC assembly: boundary DOF values are correct
   * for a constant boundary condition.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_Values_AllEqualPrescribed_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    const Real gValue = 3.14;

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(gValue));
    dbc.assemble();

    for (const auto& [local, value] : std::get<IndexMap<Real>>(dbc.getDOFs()))
      EXPECT_NEAR(value, gValue, 1e-12);
  }

  /** Complex P2 boundary values are available on every DOF owner. */
  TEST(Assembly_MPI_DirichletBC, ComplexP2ValuesReachDOFOwners)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    H1<2, Complex, Mesh<Context::MPI>> fes(std::integral_constant<size_t, 2>{}, mpiMesh);
    TrialFunction u(fes);
    const Complex prescribed(1, 2);
    DirichletBC dbc(u, ComplexFunction(prescribed));
    dbc.assemble();
    const auto& values = std::get<IndexMap<Complex>>(dbc.getDOFs());
    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);
    for (const auto& [index, value] : values)
    {
      if (begin <= index && index < end)
      {
        EXPECT_NEAR(std::abs(value - prescribed), 0, 1e-12);
      }
    }
    const size_t localCount = static_cast<size_t>(
      std::count_if(values.begin(), values.end(), [begin, end](const auto& entry) {
        return begin <= entry.first && entry.first < end;
      }));
    const size_t globalCount =
      boost::mpi::all_reduce(world, localCount, std::plus<size_t>());
    EXPECT_GT(globalCount, 0u);
  }

  // =========================================================================
  // All-geometry MPIIteration tests — Segment (1D), all supported 3D cells
  // =========================================================================

  /**
   * @brief MPIIteration over a 1D Segment mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, SegmentMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, { 10 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    // At least one cell on every rank for a 10-segment mesh
    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief Global cell count for a 1D Segment mesh equals the expected value.
   */
  TEST(Assembly_MPI_Iteration, GlobalCellCount_MatchesLocal_Segment)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, { 10 });

    // Count only owned cells on this rank (ghost cells must be excluded to
    // avoid double-counting when reducing across ranks).
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();
    size_t localCount = 0;
    {
      Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
      for (auto it = iter.getIterator(); it; ++it)
      {
        if (shard.isOwned(D, it->getIndex()))
          ++localCount;
      }
    }

    size_t totalCount = 0;
    boost::mpi::reduce(world, localCount, totalCount, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      // A UniformGrid Segment {10} has 9 cells
      EXPECT_EQ(totalCount, 9u);
    }
  }

  /**
   * @brief MPIIteration over every supported 3D cell mesh yields cells.
   */
  TEST(Assembly_MPI_Iteration, All3DMeshes_HaveCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    for (auto type : {Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(polytopeName(type));
      auto mpiMesh = distributeFromRoot(ctx, type, {4, 3, 3});

      size_t localCount = 0;
      Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
      for (auto it = iter.getIterator(); it; ++it)
        ++localCount;

      EXPECT_GT(localCount, 0u);
    }
  }

  // =========================================================================
  // All-geometry DirichletBC MPI tests
  // =========================================================================

  /**
   * @brief MPI DirichletBC assembly on Segment (1D) mesh: at least one
   * boundary DOF is found and its value matches the prescribed constant.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_Segment)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    const Real gValue = 2.71;

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, { 10 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(gValue));
    dbc.assemble();

    const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();
    size_t globalFixed = 0;
    boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      EXPECT_GT(globalFixed, 0u);
    }

    for (const auto& [local, value] : std::get<IndexMap<Real>>(dbc.getDOFs()))
      EXPECT_NEAR(value, gValue, 1e-12);
  }

  /**
   * @brief MPI DirichletBC assembly on every 3D cell fixes boundary DOFs.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_NonEmpty_All3D)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    for (auto type : {Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(polytopeName(type));

      auto mpiMesh = distributeFromRoot(ctx, type, {4, 3, 3});

      P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
      TrialFunction u(fes);

      DirichletBC dbc(u, RealFunction(1.0));
      dbc.assemble();

      const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();

      size_t globalFixed = 0;
      boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

      if (world.rank() == 0)
      {
        EXPECT_GT(globalFixed, 0u);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// main() — initializes MPI environment used by all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
