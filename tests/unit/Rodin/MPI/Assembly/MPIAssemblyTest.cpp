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
 * - P1, H1 orders 1–6, P0g, and real/complex/vector traces reach DOF owners
 * - P0, P0g, P1, and H1 orders 1–6 have non-collective fixed-layout metadata
 *   across real/complex scalar/vector ranges, including Point and empty shards
 * - Positive-dimensional cell, boundary, and nested SubMeshes preserve logical
 *   ancestry and shared DOF indices for the same families and value ranges,
 *   with exact P1 vertex-field, P0 cell-field, and P0g constant restrictions
 *
 * Run with mpirun -n 1/2/3/4/8 as registered in CMakeLists.txt.
 */
#include <algorithm>
#include <numeric>
#include <set>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

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
static boost::mpi::environment* g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /** Complex traces use the typed callable interface, not the real component pack. */
  static auto complexVectorTrace(Complex first)
  {
    return VectorFunction(size_t{3}, [first](const Point&) {
      return Math::SpatialVector<Complex>{
        first, first + Complex(1, 1), first + Complex(2, 2)};
    });
  }

  /**
   * @brief Creates a local mesh with all incidences required for sharding
   * and boundary-DOF assembly.
   */
  static Mesh<Context::Local> makeShardableMesh(
    Polytope::Type type, std::initializer_list<size_t> shape)
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
    {
      if (mesh.getShard().isOwned(mesh.getDimension(), cell->getIndex()))
        for (Index dof : fes.getDOFs(mesh.getDimension(), cell->getIndex()))
          result.insert(dof);
    }
    return result;
  }

  static const char* polytopeName(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Point:
        return "Point";
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

  class MPISpaceSizeTest : public testing::TestWithParam<Polytope::Type>
  {};

  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPISpaceSizeTest,
    testing::Values(Polytope::Type::Point, Polytope::Type::Segment,
      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) { return polytopeName(info.param); });

  /** Fixed-layout space metadata can be queried and copied on one rank only. */
  TEST_P(MPISpaceSizeTest, SpaceSizeQueriesAreNoncollective)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    const size_t dim = Polytope::Traits(GetParam()).getDimension();
    auto mesh = [&] {
      if (dim == 0)
      {
        Shard::Builder builder;
        builder.initialize(0, 3);
        if (world.rank() == world.size() - 1)
          builder.vertex(0, Math::SpatialPoint{0, 0, 0}, Shard::State::Owned);
        return Mesh<Context::MPI>::Builder(ctx).initialize(builder.finalize()).finalize();
      }
      if (dim == 1)
        return distributeFromRoot(ctx, GetParam(), {2});
      if (dim == 2)
        return distributeFromRoot(ctx, GetParam(), {2, 2});
      return distributeFromRoot(ctx, GetParam(), {2, 2, 2});
    }();
    for (size_t d = 0; d <= dim; ++d)
      for (size_t dp = 0; dp <= dim; ++dp)
        mesh.getConnectivity().compute(d, dp);
    for (size_t d = 1; d < dim; ++d)
      mesh.reconcile(d);

    const auto check = [&](const auto& fes) {
      Index begin = 0, end = 0;
      fes.getOwnershipRange(begin, end);
      const size_t expected = boost::mpi::all_reduce(
        world, static_cast<size_t>(end - begin), std::plus<size_t>());
      // Only the independent ownership-count oracle above is collective.
      // A hidden collective in any query below must fail the bounded test.
      if (world.rank() == 0)
      {
        EXPECT_EQ(fes.getSize(), expected);
        auto copied = fes;
        EXPECT_EQ(copied.getSize(), expected);
        auto moved = std::move(copied);
        EXPECT_EQ(moved.getSize(), expected);
      }
      world.barrier();
    };
    const auto range = [&]<class Range>() {
      constexpr bool scalar =
        std::is_same_v<Range, Real> || std::is_same_v<Range, Complex>;
      constexpr bool matrix = std::is_same_v<Range, Math::SpatialMatrix<Real>> ||
        std::is_same_v<Range, Math::SpatialMatrix<Complex>>;
      if constexpr (scalar)
      {
        check(P0<Range, decltype(mesh)>(mesh));
        check(P0g<Range, decltype(mesh)>(mesh));
        check(P1<Range, decltype(mesh)>(mesh));
      }
      else if constexpr (matrix)
      {
        check(P0<Range, decltype(mesh)>(mesh, 2, 3));
        check(P0g<Range, decltype(mesh)>(mesh, 2, 3));
        check(P1<Range, decltype(mesh)>(mesh, 2, 3));
      }
      else
      {
        check(P0<Range, decltype(mesh)>(mesh, 3));
        check(P0g<Range, decltype(mesh)>(mesh, 3));
        check(P1<Range, decltype(mesh)>(mesh, 3));
      }
      const auto order = [&]<size_t K>() {
        if constexpr (scalar)
          check(H1<K, Range, decltype(mesh)>(std::integral_constant<size_t, K>{}, mesh));
        else if constexpr (matrix)
          check(H1<K, Range, decltype(mesh)>(
            std::integral_constant<size_t, K>{}, mesh, 2, 3));
        else
          check(
            H1<K, Range, decltype(mesh)>(std::integral_constant<size_t, K>{}, mesh, 3));
      };
      order.template operator()<1>();
      order.template operator()<2>();
      order.template operator()<3>();
      order.template operator()<4>();
      order.template operator()<5>();
      order.template operator()<6>();
    };
    range.template operator()<Real>();
    range.template operator()<Complex>();
    range.template operator()<Math::SpatialVector<Real>>();
    range.template operator()<Math::SpatialVector<Complex>>();
    range.template operator()<Math::SpatialMatrix<Real>>();
    range.template operator()<Math::SpatialMatrix<Complex>>();
  }

  /** Reverse indices and shared-entity numbering survive unordered ghost exchange. */
  TEST_P(MPITraceGeometryTest, P0AndVectorP1GhostMapsAreBijective)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = makeMesh(ctx);

    const auto check = [&](const auto& fes, size_t entityDim) {
      const Index localCount = static_cast<Index>(fes.getShard().getSize());
      for (Index local = 0; local < localCount; ++local)
        EXPECT_EQ(fes.getLocalIndex(fes.getGlobalIndex(local)), Optional<Index>(local));

      using Entry = std::pair<Index, std::vector<Index>>;
      std::vector<Entry> localEntries;
      for (Index i = 0; i < mesh.getShard().getPolytopeCount(entityDim); ++i)
      {
        const auto& dofs = fes.getDOFs(entityDim, i);
        localEntries.emplace_back(mesh.getGlobalIndex(entityDim, i),
          std::vector<Index>(dofs.begin(), dofs.end()));
      }
      std::vector<std::vector<Entry>> gathered;
      boost::mpi::all_gather(world, localEntries, gathered);
      UnorderedMap<Index, std::vector<Index>> expected;
      for (const auto& entries : gathered)
      {
        for (const auto& [entity, dofs] : entries)
        {
          const auto [it, inserted] = expected.emplace(entity, dofs);
          if (!inserted)
          {
            EXPECT_EQ(it->second, dofs);
          }
        }
      }
    };

    P0<Real, Mesh<Context::MPI>> scalarP0(mesh);
    check(scalarP0, mesh.getDimension());
    P0<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP0(mesh, 2);
    check(vectorP0, mesh.getDimension());
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP1(mesh, 3);
    check(vectorP1, 0);
  }

  /**
   * Sparse partitions exercise empty ranks without assigning nodal meaning to
   * modal coefficients. Real/complex scalar and vector spaces must share the
   * same logical entity indices at orders one through six.
   */
  TEST_P(MPITraceGeometryTest, SparseHighOrderValueTypesShareLogicalIndices)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    const size_t dimension = Polytope::Traits(GetParam()).getDimension();
    auto mesh = dimension == 1 ? distributeFromRoot(ctx, GetParam(), {2})
      : dimension == 2         ? distributeFromRoot(ctx, GetParam(), {2, 2})
                               : distributeFromRoot(ctx, GetParam(), {2, 2, 2});
    const auto check = [&](const auto& space) {
      Index begin, end;
      space.getOwnershipRange(begin, end);
      std::vector<std::pair<Index, Index>> ranges;
      boost::mpi::all_gather(world, std::pair{begin, end}, ranges);
      Index next = 0;
      for (const auto& range : ranges)
      {
        EXPECT_EQ(range.first, next);
        EXPECT_GE(range.second, range.first);
        next = range.second;
      }
      EXPECT_EQ(next, space.getSize());
      for (Index local = 0; local < space.getShard().getSize(); ++local)
        EXPECT_EQ(
          space.getLocalIndex(space.getGlobalIndex(local)), Optional<Index>(local));
      using Record = std::pair<Index, std::vector<Index>>;
      for (size_t d = 0; d <= dimension; ++d)
      {
        std::vector<Record> local;
        for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
        {
          const auto& dofs = space.getDOFs(d, entity);
          std::vector<Index> indices(dofs.begin(), dofs.end());
          std::sort(indices.begin(), indices.end());
          local.emplace_back(mesh.getGlobalIndex(d, entity), std::move(indices));
        }
        std::vector<std::vector<Record>> gathered;
        boost::mpi::all_gather(world, local, gathered);
        IndexMap<std::vector<Index>> expected;
        for (const auto& records : gathered)
          for (const auto& [entity, indices] : records)
          {
            const auto [it, inserted] = expected.emplace(entity, indices);
            if (!inserted)
            {
              EXPECT_EQ(it->second, indices) << "dimension=" << d << " entity=" << entity;
            }
          }
      }
    };
    const auto orders = [&]<size_t K>() {
      SCOPED_TRACE(K);
      H1<K, Real, Mesh<Context::MPI>> real(std::integral_constant<size_t, K>{}, mesh);
      H1<K, Complex, Mesh<Context::MPI>> complex(
        std::integral_constant<size_t, K>{}, mesh);
      H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> realVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      EXPECT_EQ(real.getSize(), complex.getSize());
      EXPECT_EQ(realVector.getSize(), 3 * real.getSize());
      EXPECT_EQ(realVector.getSize(), complexVector.getSize());
      check(real);
      check(complex);
      check(realVector);
      check(complexVector);
      for (size_t d = 0; d <= dimension; ++d)
        for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
        {
          EXPECT_EQ(real.getDOFs(d, entity).size(), complex.getDOFs(d, entity).size());
          if (real.getDOFs(d, entity).size() == complex.getDOFs(d, entity).size())
          {
            EXPECT_TRUE((real.getDOFs(d, entity) == complex.getDOFs(d, entity)).all());
          }
          EXPECT_EQ(realVector.getDOFs(d, entity).size(),
            complexVector.getDOFs(d, entity).size());
          if (realVector.getDOFs(d, entity).size() ==
            complexVector.getDOFs(d, entity).size())
          {
            EXPECT_TRUE(
              (realVector.getDOFs(d, entity) == complexVector.getDOFs(d, entity)).all());
          }
        }
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
    orders.template operator()<5>();
    orders.template operator()<6>();
  }

  /**
   * A point is a zero-dimensional cell, not an h-refinement hierarchy.
   * Its exact constants test space construction and evaluation independently
   * of mesh coordinates. Entity ownership and global-constant DOF ownership
   * have separate contracts, including ranks without a selected local cell.
   */
  TEST_P(MPITraceGeometryTest, PointSubMeshReproducesAllValueTypes)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    const size_t dimension = Polytope::Traits(GetParam()).getDimension();
    auto parent = dimension == 1 ? distributeFromRoot(ctx, GetParam(), {5})
      : dimension == 2           ? distributeFromRoot(ctx, GetParam(), {5, 5})
                                 : distributeFromRoot(ctx, GetParam(), {5, 5, 5});
    Index localMaximum = 0;
    for (Index vertex = 0; vertex < parent.getShard().getVertexCount(); ++vertex)
      if (parent.getShard().isOwned(0, vertex))
        localMaximum = std::max(localMaximum, parent.getGlobalIndex(0, vertex));
    const Index selected = boost::mpi::all_reduce(
      world, localMaximum, [](Index a, Index b) { return std::max(a, b); });
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(parent);
    for (Index vertex = 0; vertex < parent.getShard().getVertexCount(); ++vertex)
      if (parent.getShard().isOwned(0, vertex) &&
        parent.getGlobalIndex(0, vertex) == selected)
        builder.include(0, vertex);
    auto sub = builder.finalize();
    const Mesh<Context::MPI>& mesh = sub;
    EXPECT_EQ(mesh.getDimension(), 0u);
    EXPECT_EQ(mesh.getSpaceDimension(), dimension);
    EXPECT_EQ(mesh.getPolytopeCount(0), 1u);
    size_t localOwned = 0;
    for (Index vertex = 0; vertex < mesh.getShard().getVertexCount(); ++vertex)
    {
      localOwned += mesh.getShard().isOwned(0, vertex);
      EXPECT_EQ(mesh.getGlobalIndex(0, vertex), selected);
      const Index parentLocal = sub.getPolytopeMap(0).left.at(vertex);
      EXPECT_EQ(sub.getPolytopeMap(0).right.at(parentLocal), vertex);
      EXPECT_EQ(parent.getGlobalIndex(0, parentLocal), selected);
    }
    EXPECT_EQ(boost::mpi::all_reduce(world, localOwned, std::plus<size_t>()), 1u);
    const int owner =
      boost::mpi::all_reduce(world, localOwned ? world.rank() + 1 : 0, std::plus<int>()) -
      1;
    std::vector<int> present;
    boost::mpi::all_gather(world, mesh.getShard().getVertexCount() ? 1 : 0, present);
    IndexSet holders;
    for (size_t rank = 0; rank < present.size(); ++rank)
      if (present[rank] && static_cast<int>(rank) != owner)
        holders.insert(rank);
    for (Index vertex = 0; vertex < mesh.getShard().getVertexCount(); ++vertex)
    {
      const auto& shard = mesh.getShard();
      if (shard.isOwned(0, vertex))
      {
        const auto halo = shard.getHalo(0).find(vertex);
        if (halo == shard.getHalo(0).end())
          EXPECT_TRUE(holders.empty());
        else
          EXPECT_EQ(halo->second, holders);
      }
      else
      {
        EXPECT_EQ(shard.getOwner(0).at(vertex), owner);
        EXPECT_EQ(shard.getState(0).at(vertex), Shard::State::Ghost);
      }
    }
    const auto check = [&](const auto& space) {
      using Space = std::remove_cvref_t<decltype(space)>;
      using Scalar = typename Space::ScalarType;
      SCOPED_TRACE(typeid(Space).name());
      constexpr bool vector =
        std::is_same_v<typename Space::RangeType, Math::SpatialVector<Scalar>>;
      constexpr bool matrix =
        std::is_same_v<typename Space::RangeType, Math::SpatialMatrix<Scalar>>;
      const size_t count = matrix ? 6 : vector ? 3 : 1;
      EXPECT_EQ(space.getSize(), count);
      if constexpr (matrix)
      {
        EXPECT_EQ(space.getRows(), 2);
        EXPECT_EQ(space.getColumns(), 3);
      }
      Index begin, end;
      space.getOwnershipRange(begin, end);
      std::vector<std::pair<Index, Index>> ranges;
      boost::mpi::all_gather(world, std::pair{begin, end}, ranges);
      Index next = 0;
      for (const auto& range : ranges)
      {
        EXPECT_EQ(range.first, next);
        EXPECT_GE(range.second, range.first);
        next = range.second;
      }
      EXPECT_EQ(next, count);
      const Scalar first = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Complex(2, 3);
        else
          return Real(2);
      }();
      const auto exact = [&] {
        if constexpr (matrix)
          return MatrixFunction(size_t{2}, size_t{3}, [first](const Point&) {
            Math::SpatialMatrix<Scalar> value(2, 3);
            for (size_t r = 0; r < 2; ++r)
              for (size_t s = 0; s < 3; ++s)
                value(r, s) = first + Scalar(3 * r + s);
            return value;
          });
        else if constexpr (vector)
          return VectorFunction(size_t{3}, [first](const Point&) {
            return Math::SpatialVector<Scalar>{
              first, first + Scalar(1), first + Scalar(2)};
          });
        else if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(first);
        else
          return RealFunction(first);
      }();
      GridFunction field(space);
      field = exact;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& dofs = space.getDOFs(0, cell->getIndex());
        EXPECT_EQ(dofs.size(), count);
        for (Eigen::Index component = 0; component < dofs.size(); ++component)
          EXPECT_EQ(dofs(component), component);
        const Point point(*cell, Math::SpatialPoint::Zero(0));
        const auto value = field(point);
        if constexpr (matrix)
        {
          EXPECT_EQ(value.rows(), 2);
          EXPECT_EQ(value.cols(), 3);
          if (value.rows() == 2 && value.cols() == 3)
          {
            for (size_t r = 0; r < 2; ++r)
              for (size_t s = 0; s < 3; ++s)
                EXPECT_EQ(value(r, s), first + Scalar(3 * r + s));
          }
        }
        else if constexpr (vector)
          for (size_t component = 0; component < count; ++component)
            EXPECT_EQ(value(component), first + Scalar(component));
        else
          EXPECT_EQ(value, first);
      }
      // Check provenance before restriction. Reduce the result so an empty
      // rank cannot proceed into later collectives while a holder throws.
      bool correctProvenance = true;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto pullback =
          space.getPullback({0, cell->getIndex()}, [&](const Point& point) {
            return point.getPolytope().getMesh() == static_cast<const MeshBase&>(mesh);
          });
        correctProvenance &= pullback(Math::SpatialPoint::Zero(0));
      }
      correctProvenance =
        boost::mpi::all_reduce(world, correctProvenance, std::logical_and<bool>());
      EXPECT_TRUE(correctProvenance);
      if (!correctProvenance)
        return;
      // The parent field is continuous P1: its vertex trace is unambiguous,
      // unlike a general discontinuous parent P0 trace at a shared vertex.
      const auto parentSpace = [&] {
        if constexpr (matrix)
          return P1<typename Space::RangeType, Mesh<Context::MPI>>(parent, 2, 3);
        else if constexpr (vector)
          return P1<typename Space::RangeType, Mesh<Context::MPI>>(parent, size_t{3});
        else
          return P1<typename Space::RangeType, Mesh<Context::MPI>>(parent);
      }();
      const auto scalarTrace = [first](const Point& point) {
        Scalar value = first;
        const auto& coordinates = point.getPhysicalCoordinates();
        for (size_t d = 0; d < coordinates.size(); ++d)
          value += Scalar((d + 1) * coordinates(d));
        return value;
      };
      const auto affine = [&] {
        if constexpr (matrix)
          return MatrixFunction(size_t{2}, size_t{3}, [scalarTrace](const Point& point) {
            const Scalar base = scalarTrace(point);
            Math::SpatialMatrix<Scalar> value(2, 3);
            for (size_t r = 0; r < 2; ++r)
              for (size_t s = 0; s < 3; ++s)
                value(r, s) = base + Scalar(3 * r + s);
            return value;
          });
        else if constexpr (vector)
          return VectorFunction(size_t{3}, [scalarTrace](const Point& point) {
            const Scalar base = scalarTrace(point);
            return Math::SpatialVector<Scalar>{base, base + Scalar(1), base + Scalar(2)};
          });
        else if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(scalarTrace);
        else
          return RealFunction(scalarTrace);
      }();
      GridFunction parentField(parentSpace);
      parentField = affine;
      const auto checkTrace = [&] {
        EXPECT_EQ(&field.getFiniteElementSpace(), &space);
        EXPECT_EQ(field.getData().size(), space.getSize());
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const Point point(*cell, Math::SpatialPoint::Zero(0));
          const auto actual = field(point), expected = affine(point);
          if constexpr (matrix)
          {
            EXPECT_EQ(actual.rows(), 2);
            EXPECT_EQ(actual.cols(), 3);
            if (actual.rows() == 2 && actual.cols() == 3)
            {
              for (size_t r = 0; r < 2; ++r)
                for (size_t s = 0; s < 3; ++s)
                  EXPECT_EQ(actual(r, s), expected(r, s));
            }
          }
          else if constexpr (vector)
            for (size_t component = 0; component < count; ++component)
              EXPECT_EQ(actual(component), expected(component));
          else
            EXPECT_EQ(actual, expected);
        }
      };
      // The selected entity owner has an actual point to evaluate. Other ranks
      // wait outside the local restriction/evaluation operation, so a hidden
      // collective cannot be satisfied by an all-rank call.
      if (world.rank() == owner)
      {
        field = parentField;
        checkTrace();
      }
      world.barrier();
      field = parentField;
      checkTrace();
    };
    const auto valueTypes = [&]<template <class, class> class Family>() {
      Family<Real, Mesh<Context::MPI>> real(mesh);
      Family<Complex, Mesh<Context::MPI>> complex(mesh);
      Family<Math::SpatialVector<Real>, Mesh<Context::MPI>> realVector(mesh, 3);
      Family<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(mesh, 3);
      Family<Math::SpatialMatrix<Real>, Mesh<Context::MPI>> realMatrix(mesh, 2, 3);
      Family<Math::SpatialMatrix<Complex>, Mesh<Context::MPI>> complexMatrix(mesh, 2, 3);
      check(real);
      check(complex);
      check(realVector);
      check(complexVector);
      check(realMatrix);
      check(complexMatrix);
    };
    valueTypes.template operator()<P0>();
    valueTypes.template operator()<P0g>();
    valueTypes.template operator()<P1>();
    const auto orders = [&]<size_t K>() {
      H1<K, Real, Mesh<Context::MPI>> real(std::integral_constant<size_t, K>{}, mesh);
      H1<K, Complex, Mesh<Context::MPI>> complex(
        std::integral_constant<size_t, K>{}, mesh);
      H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> realVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      H1<K, Math::SpatialMatrix<Real>, Mesh<Context::MPI>> realMatrix(
        std::integral_constant<size_t, K>{}, mesh, 2, 3);
      H1<K, Math::SpatialMatrix<Complex>, Mesh<Context::MPI>> complexMatrix(
        std::integral_constant<size_t, K>{}, mesh, 2, 3);
      check(real);
      check(complex);
      check(realVector);
      check(complexVector);
      check(realMatrix);
      check(complexMatrix);
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
    orders.template operator()<5>();
    orders.template operator()<6>();
  }

  /**
   * Positive-dimensional cell/boundary SubMeshes preserve logical ancestry
   * and shared DOF identity for all scalar/vector ranges and H1 orders 1--6.
   * Real/complex scalar/vector/matrix P1, full-dimensional P0, and P0g restrictions
   * preserve vertex/cell-label values or constants and destination layouts,
   * including nested ancestry, without tolerances. No ambiguous P0 boundary
   * trace is introduced. Values are checked on held entities, not nonexistent
   * evaluation points of empty shards.
   * Gathers below certify global ownership and agreement among holders;
   * local entity maps and round trips themselves require no collective.
   */
  TEST_P(MPITraceGeometryTest, PositiveSubMeshesPreserveAllValueTypeIndices)
  {
    const auto& world = *g_world;
    Context::MPI context(*g_env, world);
    const size_t dimension = Polytope::Traits(GetParam()).getDimension();
    auto parent = dimension == 1 ? distributeFromRoot(context, GetParam(), {2})
      : dimension == 2           ? distributeFromRoot(context, GetParam(), {2, 2})
                                 : distributeFromRoot(context, GetParam(), {2, 2, 2});
    const auto check = [&](const SubMesh<Context::MPI>& sub,
                         const std::vector<Index>& selectedParentEntities) {
      const Mesh<Context::MPI>& mesh = sub;
      const auto& immediateParent = sub.getParent();
      const size_t subDimension = mesh.getDimension();
      // Compare the explicit extraction request with the result, independently
      // of the child's ownership/DOF tables. Otherwise a dropped entity could
      // disappear from every subsequent ownership and ancestry assertion.
      std::vector<Index> actualParents;
      const auto& cellParents = sub.getPolytopeMap(subDimension).left;
      for (auto cell = mesh.getCell(); cell; ++cell)
        if (mesh.getShard().isOwned(subDimension, cell->getIndex()))
        {
          EXPECT_LT(cell->getIndex(), cellParents.size());
          if (cell->getIndex() < cellParents.size())
            actualParents.push_back(immediateParent.getGlobalIndex(
              subDimension, cellParents[cell->getIndex()]));
        }
      using Selection = std::pair<std::vector<Index>, std::vector<Index>>;
      std::vector<Selection> selections;
      boost::mpi::all_gather(
        world, Selection{selectedParentEntities, actualParents}, selections);
      std::vector<Index> expectedSelection, actualSelection;
      for (const auto& selection : selections)
      {
        expectedSelection.insert(
          expectedSelection.end(), selection.first.begin(), selection.first.end());
        actualSelection.insert(
          actualSelection.end(), selection.second.begin(), selection.second.end());
      }
      std::sort(expectedSelection.begin(), expectedSelection.end());
      std::sort(actualSelection.begin(), actualSelection.end());
      EXPECT_EQ(actualSelection, expectedSelection);
      for (size_t d = 0; d <= subDimension; ++d)
        for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
        {
          const auto& map = sub.getPolytopeMap(d);
          EXPECT_EQ(map.left.size(), mesh.getShard().getPolytopeCount(d));
          EXPECT_LT(entity, map.left.size());
          if (entity >= map.left.size())
            continue;
          const Index ancestor = map.left[entity];
          const auto inverse = map.right.find(ancestor);
          EXPECT_NE(inverse, map.right.end());
          if (inverse != map.right.end())
          {
            EXPECT_EQ(inverse->second, entity);
          }
          const auto child = mesh.getPolytope(d, entity);
          const auto parentEntity = immediateParent.getPolytope(d, ancestor);
          const auto& childVertices = child->getVertices();
          Polytope::Key mappedVertices(childVertices.size());
          bool complete = true;
          for (size_t j = 0; j < childVertices.size(); ++j)
          {
            const auto& vertices = sub.getPolytopeMap(0).left;
            EXPECT_LT(childVertices(j), vertices.size());
            if (childVertices(j) >= vertices.size())
              complete = false;
            else
              mappedVertices(j) = vertices[childVertices(j)];
          }
          if (complete)
          {
            EXPECT_TRUE(Polytope::Key::SymmetricEquality{}(
              mappedVertices, parentEntity->getVertices()));
          }
        }
      // Entity ownership is a global invariant, distinct from DOF ownership.
      for (size_t d = 0; d <= subDimension; ++d)
      {
        std::vector<std::pair<Index, Index>> ancestry;
        std::vector<Index> owned;
        for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
        {
          const Index global = mesh.getGlobalIndex(d, entity);
          const auto& ancestors = sub.getPolytopeMap(d).left;
          EXPECT_LT(entity, ancestors.size());
          if (entity < ancestors.size())
            ancestry.emplace_back(
              global, immediateParent.getGlobalIndex(d, ancestors[entity]));
          if (mesh.getShard().isOwned(d, entity))
            owned.push_back(global);
        }
        using EntityState =
          std::pair<std::vector<std::pair<Index, Index>>, std::vector<Index>>;
        std::vector<EntityState> gathered;
        boost::mpi::all_gather(world, EntityState{ancestry, owned}, gathered);
        IndexMap<Index> expectedAncestor;
        IndexMap<int> owners;
        IndexMap<IndexSet> holders;
        for (size_t rank = 0; rank < gathered.size(); ++rank)
        {
          for (const auto& [entity, ancestor] : gathered[rank].first)
          {
            const auto [it, inserted] = expectedAncestor.emplace(entity, ancestor);
            if (!inserted)
            {
              EXPECT_EQ(it->second, ancestor);
            }
            holders[entity].insert(rank);
          }
          for (const Index entity : gathered[rank].second)
          {
            const auto [it, inserted] = owners.emplace(entity, static_cast<int>(rank));
            EXPECT_TRUE(inserted) << "dimension=" << d << " entity=" << entity;
          }
        }
        for (const auto& [entity, ancestor] : expectedAncestor)
          EXPECT_TRUE(owners.contains(entity))
            << "dimension=" << d << " entity=" << entity;
        for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
        {
          const Index global = mesh.getGlobalIndex(d, entity);
          const auto owner = owners.find(global);
          if (owner == owners.end())
            continue;
          const auto& shard = mesh.getShard();
          if (shard.isOwned(d, entity))
          {
            EXPECT_EQ(owner->second, world.rank());
            auto expectedHalo = holders[global];
            expectedHalo.erase(world.rank());
            const auto halo = shard.getHalo(d).find(entity);
            if (halo == shard.getHalo(d).end())
              EXPECT_TRUE(expectedHalo.empty());
            else
              EXPECT_EQ(halo->second, expectedHalo);
          }
          else
          {
            const auto ghostOwner = shard.getOwner(d).find(entity);
            EXPECT_NE(ghostOwner, shard.getOwner(d).end());
            if (ghostOwner != shard.getOwner(d).end())
            {
              EXPECT_EQ(ghostOwner->second, owner->second);
            }
          }
        }
      }
      const auto checkSpace = [&](const auto& space) {
        using Space = std::remove_cvref_t<decltype(space)>;
        constexpr bool cellOnly =
          std::is_same_v<Space, P0<typename Space::RangeType, Mesh<Context::MPI>>>;
        SCOPED_TRACE(typeid(space).name());
        // Exercise the point provenance used by parent GridFunction evaluation.
        // This is strictly rank-local ancestry, not a physical-coordinate match.
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const Index child = cell->getIndex();
          const auto& ancestors = sub.getPolytopeMap(subDimension).left;
          EXPECT_LT(child, ancestors.size());
          if (child >= ancestors.size())
            continue; // The entity-map checks above report the missing entry.
          const Index ancestor = ancestors[child];
          Index rootAncestor = ancestor;
          if (immediateParent.isSubMesh())
          {
            const auto& parentMap =
              immediateParent.asSubMesh().getPolytopeMap(subDimension);
            EXPECT_LT(rootAncestor, parentMap.left.size());
            if (rootAncestor >= parentMap.left.size())
              continue;
            rootAncestor = parentMap.left[rootAncestor];
          }
          const auto pullback =
            space.getPullback({subDimension, child}, [&](const Point& point) {
              if (point.getPolytope().getMesh() != static_cast<const MeshBase&>(mesh))
                return false;
              const auto immediate = immediateParent.inclusion(point);
              const auto root = parent.inclusion(point);
              return immediate && root &&
                immediate->getPolytope().getMesh() ==
                static_cast<const MeshBase&>(immediateParent) &&
                immediate->getPolytope().getDimension() == subDimension &&
                immediate->getPolytope().getIndex() == ancestor &&
                root->getPolytope().getMesh() == static_cast<const MeshBase&>(parent) &&
                root->getPolytope().getDimension() == subDimension &&
                root->getPolytope().getIndex() == rootAncestor;
            });
          EXPECT_TRUE(pullback(Math::SpatialPoint::Zero(subDimension)));
        }
        for (Index local = 0; local < space.getShard().getSize(); ++local)
        {
          if constexpr (requires { space.getLocalIndex(Index{}); })
            EXPECT_EQ(
              space.getLocalIndex(space.getGlobalIndex(local)), Optional<Index>(local));
          else
            // P0g has a replicated constant view, not a local inverse-map API.
            EXPECT_EQ(space.getGlobalIndex(local), local);
        }
        Index begin, end;
        space.getOwnershipRange(begin, end);
        std::vector<std::pair<Index, Index>> ranges;
        boost::mpi::all_gather(world, std::pair{begin, end}, ranges);
        Index next = 0;
        for (const auto& range : ranges)
        {
          EXPECT_EQ(range.first, next);
          EXPECT_GE(range.second, range.first);
          next = range.second;
        }
        EXPECT_EQ(next, space.getSize());
        using Record = std::pair<Index, std::vector<Index>>;
        // P0 DOFs are defined on cells, not on their lower-dimensional faces.
        for (size_t d = cellOnly ? subDimension : 0; d <= subDimension; ++d)
        {
          std::vector<Record> local;
          for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
          {
            const auto& dofs = space.getDOFs(d, entity);
            std::vector<Index> indices(dofs.begin(), dofs.end());
            std::sort(indices.begin(), indices.end());
            for (const Index dof : indices)
            {
              EXPECT_LT(dof, space.getSize());
              if constexpr (requires { space.getLocalIndex(Index{}); })
                EXPECT_TRUE(space.getLocalIndex(dof).has_value());
              else
                EXPECT_LT(dof, space.getShard().getSize());
            }
            local.emplace_back(mesh.getGlobalIndex(d, entity), std::move(indices));
          }
          std::vector<std::vector<Record>> gathered;
          boost::mpi::all_gather(world, local, gathered);
          IndexMap<std::vector<Index>> expected;
          for (const auto& records : gathered)
            for (const auto& [entity, indices] : records)
            {
              const auto [it, inserted] = expected.emplace(entity, indices);
              if (!inserted)
              {
                EXPECT_EQ(it->second, indices)
                  << "dimension=" << d << " entity=" << entity;
              }
            }
        }
      };
      const auto compareRanges = [&](const auto& real, const auto& complex,
                                   const auto& realVector, const auto& complexVector) {
        constexpr bool cellOnly = std::is_same_v<std::remove_cvref_t<decltype(real)>,
          P0<Real, Mesh<Context::MPI>>>;
        EXPECT_EQ(real.getSize(), complex.getSize());
        EXPECT_EQ(realVector.getSize(), 3 * real.getSize());
        EXPECT_EQ(realVector.getSize(), complexVector.getSize());
        checkSpace(real);
        checkSpace(complex);
        checkSpace(realVector);
        checkSpace(complexVector);
        for (size_t d = cellOnly ? subDimension : 0; d <= subDimension; ++d)
          for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
          {
            // getDOFs may return a reusable scratch view. Compare snapshots,
            // not two references that could alias the next query's storage.
            const IndexArray scalar = real.getDOFs(d, entity);
            const IndexArray complexScalar = complex.getDOFs(d, entity);
            EXPECT_EQ(scalar.size(), complexScalar.size());
            if (scalar.size() == complexScalar.size())
            {
              EXPECT_TRUE((scalar == complexScalar).all());
            }
            const IndexArray vector = realVector.getDOFs(d, entity);
            const IndexArray complexValues = complexVector.getDOFs(d, entity);
            EXPECT_EQ(vector.size(), complexValues.size());
            if (vector.size() == complexValues.size())
            {
              EXPECT_TRUE((vector == complexValues).all());
            }
          }
      };
      const auto compareMatrices = [&](const auto& scalar, const auto& realMatrix,
                                     const auto& complexMatrix) {
        constexpr bool cellOnly = std::is_same_v<std::remove_cvref_t<decltype(scalar)>,
          P0<Real, Mesh<Context::MPI>>>;
        EXPECT_EQ(realMatrix.getRows(), 2);
        EXPECT_EQ(realMatrix.getColumns(), 3);
        EXPECT_EQ(complexMatrix.getRows(), 2);
        EXPECT_EQ(complexMatrix.getColumns(), 3);
        EXPECT_EQ(realMatrix.getSize(), 6 * scalar.getSize());
        EXPECT_EQ(complexMatrix.getSize(), realMatrix.getSize());
        checkSpace(realMatrix);
        checkSpace(complexMatrix);
        for (size_t d = cellOnly ? subDimension : 0; d <= subDimension; ++d)
          for (Index entity = 0; entity < mesh.getShard().getPolytopeCount(d); ++entity)
          {
            const IndexArray scalarDOFs = scalar.getDOFs(d, entity);
            const IndexArray matrixDOFs = realMatrix.getDOFs(d, entity);
            const IndexArray complexDOFs = complexMatrix.getDOFs(d, entity);
            ASSERT_EQ(matrixDOFs.size(), 6 * scalarDOFs.size());
            ASSERT_EQ(complexDOFs.size(), matrixDOFs.size());
            for (Eigen::Index a = 0; a < matrixDOFs.size(); ++a)
            {
              EXPECT_EQ(matrixDOFs(a), 6 * scalarDOFs(a / 6) + a % 6);
              EXPECT_EQ(complexDOFs(a), matrixDOFs(a));
            }
          }
      };
      const auto families = [&]<template <class, class> class Family>() {
        Family<Real, Mesh<Context::MPI>> real(mesh);
        Family<Complex, Mesh<Context::MPI>> complex(mesh);
        Family<Math::SpatialVector<Real>, Mesh<Context::MPI>> realVector(mesh, 3);
        Family<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(mesh, 3);
        compareRanges(real, complex, realVector, complexVector);
        Family<Math::SpatialMatrix<Real>, Mesh<Context::MPI>> realMatrix(mesh, 2, 3);
        Family<Math::SpatialMatrix<Complex>, Mesh<Context::MPI>> complexMatrix(
          mesh, 2, 3);
        compareMatrices(real, realMatrix, complexMatrix);
      };
      families.template operator()<P0>();
      families.template operator()<P0g>();
      families.template operator()<P1>();
      const auto orders = [&]<size_t K>() {
        SCOPED_TRACE(K);
        H1<K, Real, Mesh<Context::MPI>> real(std::integral_constant<size_t, K>{}, mesh);
        H1<K, Complex, Mesh<Context::MPI>> complex(
          std::integral_constant<size_t, K>{}, mesh);
        H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> realVector(
          std::integral_constant<size_t, K>{}, mesh, 3);
        H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
          std::integral_constant<size_t, K>{}, mesh, 3);
        compareRanges(real, complex, realVector, complexVector);
      };
      orders.template operator()<1>();
      orders.template operator()<2>();
      orders.template operator()<3>();
      orders.template operator()<4>();
      orders.template operator()<5>();
      orders.template operator()<6>();
      const auto checkRestriction =
        [&]<class Range, template <class, class> class Family = P1>() {
          using Space = Family<Range, Mesh<Context::MPI>>;
          constexpr bool CellConstant =
            std::is_same_v<Space, P0<Range, Mesh<Context::MPI>>>;
          constexpr bool GlobalConstant =
            std::is_same_v<Space, P0g<Range, Mesh<Context::MPI>>>;
          using Scalar = typename Space::ScalarType;
          constexpr bool vector = std::is_same_v<Range, Math::SpatialVector<Scalar>>;
          constexpr bool matrix = std::is_same_v<Range, Math::SpatialMatrix<Scalar>>;
          SCOPED_TRACE(typeid(Space).name());
          const auto sourceSpace = [&] {
            if constexpr (matrix)
              return Space(parent, 2, 3);
            else if constexpr (vector)
              return Space(parent, 3);
            else
              return Space(parent);
          }();
          const auto targetSpace = [&] {
            if constexpr (matrix)
              return Space(mesh, 2, 3);
            else if constexpr (vector)
              return Space(mesh, 3);
            else
              return Space(mesh);
          }();
          GridFunction source(sourceSpace), restricted(targetSpace);
          const auto value = [](Index entity, size_t component) {
            const Scalar base = [&] {
              if constexpr (std::is_same_v<Scalar, Complex>)
                return Complex(Real(entity + 1), Real(2 * entity + 1));
              else
                return Real(entity + 1);
            }();
            return base + Scalar(component);
          };
          // P1 uses vertex labels; P0 uses cell labels and is restricted only
          // to full-dimensional cells. No incident-cell choice defines a P0
          // boundary trace here. P0g uses the same constant label on every held
          // entity. All holder coefficients are available locally.
          const size_t entityDimension = CellConstant ? dimension : 0;
          for (Index entity = 0;
               entity < parent.getShard().getPolytopeCount(entityDimension); ++entity)
          {
            const IndexArray dofs = sourceSpace.getDOFs(entityDimension, entity);
            const Index global =
              GlobalConstant ? Index{0} : parent.getGlobalIndex(entityDimension, entity);
            for (size_t component = 0; component < static_cast<size_t>(dofs.size());
                 ++component)
              source[dofs(component)] = value(global, component);
          }
          // Different spaces must interpolate, rather than copy/rebind the
          // source layout. The parent inclusion follows the actual value path.
          // First exercise native coefficient restriction on one rank only:
          // this operation has no global result and must not communicate.
          // The barrier is test-protocol synchronization outside the operation;
          // the bounded MPI test detects an accidental hidden collective.
          const auto checkValues = [&] {
            EXPECT_EQ(&restricted.getFiniteElementSpace(), &targetSpace);
            EXPECT_EQ(restricted.getData().size(), targetSpace.getSize());
            for (Index entity = 0;
                 entity < mesh.getShard().getPolytopeCount(entityDimension); ++entity)
            {
              const auto& ancestors = sub.getPolytopeMap(entityDimension).left;
              EXPECT_LT(entity, ancestors.size());
              if (entity >= ancestors.size())
                continue;
              Index ancestor = ancestors[entity];
              if (immediateParent.isSubMesh())
              {
                const auto& map =
                  immediateParent.asSubMesh().getPolytopeMap(entityDimension);
                EXPECT_LT(ancestor, map.left.size());
                if (ancestor >= map.left.size())
                  continue;
                ancestor = map.left[ancestor];
              }
              const Index global = GlobalConstant
                ? Index{0}
                : parent.getGlobalIndex(entityDimension, ancestor);
              const IndexArray dofs = targetSpace.getDOFs(entityDimension, entity);
              for (size_t component = 0; component < static_cast<size_t>(dofs.size());
                   ++component)
                EXPECT_EQ(restricted[dofs(component)], value(global, component));
            }
          };
          if (world.rank() == 0)
          {
            restricted = source;
            checkValues();
          }
          world.barrier();
          restricted = source;
          checkValues();
        };
      checkRestriction.template operator()<Real>();
      checkRestriction.template operator()<Complex>();
      checkRestriction.template operator()<Math::SpatialVector<Real>>();
      checkRestriction.template operator()<Math::SpatialVector<Complex>>();
      checkRestriction.template operator()<Math::SpatialMatrix<Real>>();
      checkRestriction.template operator()<Math::SpatialMatrix<Complex>>();
      if (subDimension == dimension)
      {
        checkRestriction.template operator()<Real, P0>();
        checkRestriction.template operator()<Complex, P0>();
        checkRestriction.template operator()<Math::SpatialVector<Real>, P0>();
        checkRestriction.template operator()<Math::SpatialVector<Complex>, P0>();
        checkRestriction.template operator()<Math::SpatialMatrix<Real>, P0>();
        checkRestriction.template operator()<Math::SpatialMatrix<Complex>, P0>();
      }
      checkRestriction.template operator()<Real, P0g>();
      checkRestriction.template operator()<Complex, P0g>();
      checkRestriction.template operator()<Math::SpatialVector<Real>, P0g>();
      checkRestriction.template operator()<Math::SpatialVector<Complex>, P0g>();
      checkRestriction.template operator()<Math::SpatialMatrix<Real>, P0g>();
      checkRestriction.template operator()<Math::SpatialMatrix<Complex>, P0g>();
    };
    for (const auto& [selectedDimension, sparse] : {std::pair{dimension, false},
           std::pair{dimension - 1, false}, std::pair{dimension, true}})
    {
      if (selectedDimension == 0)
        continue; // Point submeshes have their own exact-value coverage.
      SCOPED_TRACE(::testing::Message()
        << "submesh dimension=" << selectedDimension << " sparse=" << sparse);
      SubMesh<Context::MPI>::Builder builder;
      builder.initialize(parent);
      std::vector<Index> selected;
      if (selectedDimension == dimension)
      {
        for (auto cell = parent.getCell(); cell; ++cell)
          if (parent.getShard().isOwned(dimension, cell->getIndex()) &&
            (!sparse || parent.getGlobalIndex(dimension, cell->getIndex()) % 2 == 0))
          {
            builder.include(dimension, cell->getIndex());
            selected.push_back(parent.getGlobalIndex(dimension, cell->getIndex()));
          }
      }
      else
      {
        for (auto face = parent.getBoundary(); face; ++face)
          if (parent.getShard().isOwned(selectedDimension, face->getIndex()))
          {
            builder.include(selectedDimension, face->getIndex());
            selected.push_back(
              parent.getGlobalIndex(selectedDimension, face->getIndex()));
          }
      }
      auto sub = builder.finalize();
      EXPECT_EQ(sub.getDimension(), selectedDimension);
      check(sub, selected);
      SubMesh<Context::MPI>::Builder nestedBuilder;
      nestedBuilder.initialize(sub);
      std::vector<Index> nestedSelected;
      for (auto cell = sub.getCell(); cell; ++cell)
        if (sub.getShard().isOwned(selectedDimension, cell->getIndex()))
        {
          nestedBuilder.include(selectedDimension, cell->getIndex());
          nestedSelected.push_back(
            sub.getGlobalIndex(selectedDimension, cell->getIndex()));
        }
      auto nested = nestedBuilder.finalize();
      EXPECT_EQ(nested.getDimension(), selectedDimension);
      check(nested, nestedSelected);
    }
  }

  /** Every boundary DOF reaches its owner, independently of boundary-face ownership. */
  TEST_P(MPITraceGeometryTest, BoundaryConstraintsReachOwnersAcrossSpaces)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = makeMesh(ctx);
    auto probe = [&](const auto& fes, const auto& value, const char* name) {
      if constexpr (requires { fes.getLocalIndex(Index{}); })
      {
        const Index localSize = fes.getShard().getSize();
        Index mismatch = localSize;
        for (Index local = 0; local < localSize; ++local)
        {
          if (fes.getLocalIndex(fes.getGlobalIndex(local)) != Optional<Index>(local))
          {
            mismatch = local;
            break;
          }
        }
        EXPECT_EQ(mismatch, localSize) << name << " rank=" << world.rank();
      }
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
    P0g<Real, Mesh<Context::MPI>> p0g(mesh);
    probe(p0g, realValue, "P0g");
    P0g<Complex, Mesh<Context::MPI>> complexP0g(mesh);
    probe(complexP0g, complexValue, "P0g complex");
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP1(mesh, 3);
    probe(vectorP1, vectorValue, "P1 vector");
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP0g(mesh, 3);
    probe(vectorP0g, vectorValue, "P0g vector");
    P1<Complex, Mesh<Context::MPI>> complexP1(mesh);
    probe(complexP1, complexValue, "P1 complex");
    const auto complexVectorValue = complexVectorTrace(Complex(1, 2));
    P1<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVectorP1(mesh, 3);
    probe(complexVectorP1, complexVectorValue, "P1 complex vector");
    P0g<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVectorP0g(mesh, 3);
    probe(complexVectorP0g, complexVectorValue, "P0g complex vector");
    // Each order's metadata is released before constructing the next order.
    const auto orders = [&]<size_t K>() {
      SCOPED_TRACE(K);
      H1<K, Real, Mesh<Context::MPI>> scalar(std::integral_constant<size_t, K>{}, mesh);
      H1<K, Complex, Mesh<Context::MPI>> complex(
        std::integral_constant<size_t, K>{}, mesh);
      H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> vector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      probe(scalar, realValue, "H1 real");
      probe(complex, complexValue, "H1 complex");
      probe(vector, vectorValue, "H1 vector");
      probe(complexVector, complexVectorValue, "H1 complex vector");
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
    orders.template operator()<5>();
    orders.template operator()<6>();
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
    {
      if (parent.getShard().isOwned(D, cell->getIndex()))
        builder.include(D, cell->getIndex());
    }
    auto sub = builder.finalize();
    const Mesh<Context::MPI>& mesh = sub;
    const auto probe = [&](const auto& fes) {
      using Space = std::remove_cvref_t<decltype(fes)>;
      using Scalar = typename Space::ScalarType;
      const auto prescribed = [] {
        if constexpr (std::is_same_v<typename Space::RangeType,
                        Math::SpatialVector<Real>>)
          return VectorFunction{RealFunction(1), RealFunction(2), RealFunction(3)};
        else if constexpr (std::is_same_v<typename Space::RangeType,
                             Math::SpatialVector<Complex>>)
          return complexVectorTrace(Complex(1, 2));
        else if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(Complex(1, 2));
        else
          return RealFunction(1);
      }();
      TrialFunction u(fes);
      auto dbc = DirichletBC(u, prescribed);
      dbc.assemble();
      const auto& values = std::get<IndexMap<Scalar>>(dbc.getDOFs());
      const auto required = requiredDOFs(fes);
      std::set<Index> expected;
      for (auto face = mesh.getFace(); face; ++face)
      {
        if (physical.contains(
              mesh.getShard().getPolytopeMap(D - 1).left.at(face->getIndex())))
          for (Index dof : fes.getDOFs(D - 1, face->getIndex()))
          {
            if (required.contains(dof))
              expected.insert(dof);
          }
      }
      EXPECT_EQ(values.size(), expected.size());
      for (Index dof : expected)
        EXPECT_TRUE(values.contains(dof));
    };
    P1<Real, Mesh<Context::MPI>> p1(mesh);
    H1<2, Real, Mesh<Context::MPI>> p2(std::integral_constant<size_t, 2>{}, mesh);
    probe(p1);
    probe(p2);
    P1<Complex, Mesh<Context::MPI>> complexP1(mesh);
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP1(mesh, 3);
    P1<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVectorP1(mesh, 3);
    probe(complexP1);
    probe(vectorP1);
    probe(complexVectorP1);
    const auto orders = [&]<size_t K>() {
      H1<K, Real, Mesh<Context::MPI>> scalar(std::integral_constant<size_t, K>{}, mesh);
      H1<K, Complex, Mesh<Context::MPI>> complex(
        std::integral_constant<size_t, K>{}, mesh);
      H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> vector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      probe(scalar);
      probe(complex);
      probe(vector);
      probe(complexVector);
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
    orders.template operator()<5>();
    orders.template operator()<6>();
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
        else if constexpr (std::is_same_v<typename Space::RangeType,
                             Math::SpatialVector<Complex>>)
          return complexVectorTrace(Complex(2, 3));
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
      {
        for (Index dof : fes.getDOFs(mesh.getDimension() - 1, it->getIndex()))
          local.push_back(dof);
      }
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
    P0g<Real, Mesh<Context::MPI>> p0g(mesh);
    probe(p0g, "P0g real");
    P0g<Complex, Mesh<Context::MPI>> complexP0g(mesh);
    probe(complexP0g, "P0g complex");
    P1<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP1(mesh, 3);
    probe(vectorP1, "P1 vector");
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> vectorP0g(mesh, 3);
    probe(vectorP0g, "P0g vector");
    P1<Complex, Mesh<Context::MPI>> complexP1(mesh);
    probe(complexP1, "P1 complex");
    P1<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVectorP1(mesh, 3);
    probe(complexVectorP1, "P1 complex vector");
    P0g<Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVectorP0g(mesh, 3);
    probe(complexVectorP0g, "P0g complex vector");
    // Retain one order's four range layouts, not every preceding order.
    const auto orders = [&]<size_t K>() {
      SCOPED_TRACE(K);
      H1<K, Real, Mesh<Context::MPI>> scalar(std::integral_constant<size_t, K>{}, mesh);
      H1<K, Complex, Mesh<Context::MPI>> complex(
        std::integral_constant<size_t, K>{}, mesh);
      H1<K, Math::SpatialVector<Real>, Mesh<Context::MPI>> vector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      probe(scalar, "H1 real");
      probe(complex, "H1 complex");
      probe(vector, "H1 vector");
      probe(complexVector, "H1 complex vector");
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
    orders.template operator()<5>();
    orders.template operator()<6>();
  }

  /**
   * Matrix traces use the same required-DOF contract as scalar/vector traces.
   * The oracle classifies faces by their original logical parent IDs, not by
   * the boundary iterator used by assembly or by coordinate comparisons.
   * P0g is globally supported, including on empty shards: its prescribed
   * payload is collectively selected once and retained on every holder.
   * Other spaces use the owned-DOF/owned-cell required set. The index oracle
   * itself introduces no exchange beyond the initial parent-ID broadcast.
   */
  TEST_P(MPITraceGeometryTest, MatrixBoundaryAndIdentificationMatchLogicalTrace)
  {
    Context::MPI ctx(*g_env, *g_world);
    std::vector<Index> parentBoundary;
    const size_t D = Polytope::Traits(GetParam()).getDimension();
    // A bounded grid exercises cell/face/edge/vertex DOFs through degree six
    // without retaining the large legacy tetrahedral ownership workload.
    auto mesh = D == 1 ? distributeFromRoot(ctx, GetParam(), {3}, &parentBoundary)
      : D == 2         ? distributeFromRoot(ctx, GetParam(), {3, 3}, &parentBoundary)
                       : distributeFromRoot(ctx, GetParam(), {3, 3, 3}, &parentBoundary);
    const std::set<Index> physical(parentBoundary.begin(), parentBoundary.end());
    const auto probe = [&](const auto& fes) {
      using Space = std::remove_cvref_t<decltype(fes)>;
      using Scalar = typename Space::ScalarType;
      SCOPED_TRACE(typeid(Space).name());
      const auto required = requiredDOFs(fes);
      std::set<Index> expected;
      for (auto face = mesh.getFace(); face; ++face)
        if (physical.contains(
              mesh.getShard().getPolytopeMap(D - 1).left.at(face->getIndex())))
          for (Index dof : fes.getDOFs(D - 1, face->getIndex()))
            if (required.contains(dof))
              expected.insert(dof);
      if constexpr (std::is_same_v<typename Space::ElementType,
                      P0gElement<typename Space::RangeType>>)
      {
        // The global constant basis is not attached to a local entity.
        // Even an empty shard holds its components and receives global data.
        ASSERT_FALSE(physical.empty());
        expected.clear();
        for (Index dof = 0; dof < fes.getSize(); ++dof)
          expected.insert(dof);
      }
      Math::SpatialMatrix<Scalar> prescribed(2, 3);
      for (size_t r = 0; r < 2; ++r)
        for (size_t c = 0; c < 3; ++c)
        {
          prescribed(r, c) = Scalar(1 + 3 * r + c);
          if constexpr (std::is_same_v<Scalar, Complex>)
            prescribed(r, c) += Complex(0, 1 + r + c);
        }
      TrialFunction u(fes);
      TrialFunction v(fes);
      auto value = DirichletBC(u, MatrixFunction(prescribed));
      auto affine = DirichletBC(u, -v, MatrixFunction(prescribed));
      if constexpr (!std::is_same_v<typename Space::ElementType,
                      P0gElement<typename Space::RangeType>>)
      {
        // Construct fields collectively above, then assemble on one rank
        // only. The barrier is a test-protocol rendezvous, not reconciliation.
        if (g_world->rank() == 0)
        {
          value.assemble();
          affine.assemble();
        }
        g_world->barrier();
      }
      value.assemble();
      affine.assemble();
      const auto& values = std::get<IndexMap<Scalar>>(value.getDOFs());
      const auto& rows =
        std::get<typename DirichletBCBase<Scalar>::IdentifiedDOFs>(affine.getDOFs());
      const auto& offsets = affine.getIdentificationValues();
      EXPECT_EQ(values.size(), expected.size());
      EXPECT_EQ(rows.size(), expected.size());
      EXPECT_EQ(offsets.size(), expected.size());
      for (Index dof : expected)
      {
        EXPECT_TRUE(values.contains(dof));
        EXPECT_TRUE(offsets.contains(dof));
        const auto row = rows.find(dof);
        ASSERT_NE(row, rows.end());
        EXPECT_EQ(row->second.first.size(), row->second.second.size());
        bool matchingMaster = false;
        for (Index i = 0; i < static_cast<Index>(row->second.first.size()); ++i)
          matchingMaster |= row->second.first[i] == dof;
        EXPECT_TRUE(matchingMaster);
      }
    };
    const auto ranges = [&]<class Scalar>(std::type_identity<Scalar>) {
      using Range = Math::SpatialMatrix<Scalar>;
      P0g<Range, Mesh<Context::MPI>> p0g(mesh, 2, 3);
      probe(p0g);
      P1<Range, Mesh<Context::MPI>> p1(mesh, 2, 3);
      probe(p1);
      const auto order = [&]<size_t K>(std::integral_constant<size_t, K>) {
        H1<K, Range, Mesh<Context::MPI>> fes(
          std::integral_constant<size_t, K>{}, mesh, 2, 3);
        probe(fes);
      };
      order(std::integral_constant<size_t, 1>{});
      order(std::integral_constant<size_t, 2>{});
      order(std::integral_constant<size_t, 3>{});
      order(std::integral_constant<size_t, 4>{});
      order(std::integral_constant<size_t, 5>{});
      order(std::integral_constant<size_t, 6>{});
    };
    ranges(std::type_identity<Real>{});
    ranges(std::type_identity<Complex>{});
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
      {
        if (shard.isOwned(dim, cell->getIndex()))
          for (Index dof : space.getDOFs(dim, cell->getIndex()))
            required.insert(dof);
      }

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
      {
        for (const auto& [dof, face] : records)
          expected[dof].insert(face);
      }
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
      H1<K, Math::SpatialVector<Complex>, Mesh<Context::MPI>> complexVector(
        std::integral_constant<size_t, K>{}, mesh, 3);
      probe(complexVector, "H1 complex vector");
    };
    orders.template operator()<1>();
    orders.template operator()<2>();
    orders.template operator()<3>();
    orders.template operator()<4>();
    orders.template operator()<5>();
    orders.template operator()<6>();
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
      {
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
      }
      for (const auto& records : gathered)
      {
        for (const auto& [gid, declaration] : records)
        {
          const auto owner = owners.find(gid);
          EXPECT_NE(owner, owners.end());
          if (owner != owners.end())
          {
            EXPECT_EQ(owner->second, declaration.first);
          }
        }
      }
      for (Index i = 0; i < shard.getPolytopeCount(dim); ++i)
      {
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
    {
      if (it->getAttribute() == selected)
        ++localSelected;
    }
    size_t visibleSelected = 0;
    for (auto it = mpiMesh.getFace(); it; ++it)
    {
      if (it->getAttribute() == selected)
        ++visibleSelected;
    }
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
      {
        if (!parent.isBoundary(face->getIndex()))
          parent.setAttribute({1, face->getIndex()}, selected);
      }
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
      {
        if (mesh.getShard().isOwned(1, face->getIndex()) &&
          face->getAttribute() == selected)
          for (Index dof : fes.getDOFs(1, face->getIndex()))
            local.push_back(dof);
      }
      std::vector<std::vector<Index>> gathered;
      boost::mpi::all_gather(*g_world, local, gathered);
      const auto required = requiredDOFs(fes);
      std::set<Index> expected;
      for (const auto& indices : gathered)
      {
        for (Index dof : indices)
        {
          if (required.contains(dof))
            expected.insert(dof);
        }
      }
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
    {
      for (Index dof : slave.getDOFs(1, it->getIndex()))
        boundary.push_back(dof);
    }
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
    {
      for (Index dof : indices)
      {
        if (local.contains(dof))
          expected.insert(dof);
      }
    }
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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, {4, 4});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    P1<Real, Mesh<Context::MPI>> mpiFes(mpiMesh);
    TrialFunction uMPI(mpiFes);
    DirichletBC dbcMPI(uMPI, RealFunction(1.0));
    dbcMPI.assemble();
    const size_t mpiFixed = std::get<IndexMap<Real>>(dbcMPI.getDOFs()).size();

    // Sequential reference
    auto localMesh = makeShardableMesh(Polytope::Type::Triangle, {4, 4});
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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, {10});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, {10});

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
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, {10});

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
  g_env = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
