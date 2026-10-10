/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <Rodin/Alert.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/Variational.h>
#include <Rodin/MPI/Variational/Interpolation.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  boost::mpi::environment* g_environment = nullptr;
  boost::mpi::communicator* g_world = nullptr;

  constexpr std::array geometries = {Polytope::Type::Segment, Polytope::Type::Triangle,
    Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
    Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge};

  /** Leaves the P0g DOF owner empty whenever more than one rank is present. */
  class EmptyRootPartitioner final : public Partitioner
  {
    public:
      EmptyRootPartitioner(const Mesh<Context::Local>& mesh, size_t count)
        : m_base(mesh),
          m_count(count)
      {
        partition(count, mesh.getDimension());
      }

      const Mesh<Context::Local>& getMesh() const override
      {
        return m_base.getMesh();
      }

      void partition(size_t count, size_t dimension) override
      {
        m_count = count;
        m_base.partition(count > 1 ? count - 1 : count, dimension);
      }

      size_t getPartition(Index i) const override
      {
        return m_base.getPartition(i) + (m_count > 1);
      }

      size_t getCount() const override
      {
        return m_count;
      }

    private:
      BalancedCompactPartitioner m_base;
      size_t m_count;
  };

  Mesh<Context::MPI> distribute(const Context::MPI& context, Polytope::Type geometry)
  {
    const auto& comm = context.getCommunicator();
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    Sharder<Context::MPI> sharder(context);
    if (comm.rank() == 0)
    {
      auto mesh = dimension == 1 ? Mesh<Context::Local>::UniformGrid(geometry, {3})
        : dimension == 2         ? Mesh<Context::Local>::UniformGrid(geometry, {3, 3})
                                 : Mesh<Context::Local>::UniformGrid(geometry, {3, 3, 3});
      auto& connectivity = mesh.getConnectivity();
      connectivity.compute(dimension, dimension);
      for (size_t d = 1; d <= dimension; ++d)
      {
        connectivity.compute(dimension, d - 1);
        connectivity.compute(d, d - 1);
      }
      connectivity.compute(dimension - 1, dimension);
      EmptyRootPartitioner partitioner(mesh, comm.size());
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    return sharder.gather(0);
  }

  /**
   * Independent global oracle: every owned eligible entity contributes its
   * DOFs and functional values. The smallest distributed entity is selected
   * after gathering those records, without consulting interpolation sources.
   * All acceptance checks use exact indices and exact scalar equality.
   */
  template <class FES>
  void check(const FES& fes, Geometry::Region region, int filter)
  {
    using Scalar = typename FES::ScalarType;
    using Range = typename FES::RangeType;
    const auto& mesh = fes.getMesh();
    const auto& shard = mesh.getShard();
    const auto& comm = mesh.getContext().getCommunicator();
    SCOPED_TRACE(::testing::Message()
      << "rank=" << comm.rank() << " region=" << static_cast<int>(region)
      << " filter=" << filter);
    const size_t dimension = mesh.getDimension();
    const size_t sourceDimension =
      region == Geometry::Region::Cells ? dimension : dimension - 1;
    const auto predicate = [&](const Polytope& entity) {
      const Index source = mesh.getGlobalIndex(entity.getDimension(), entity.getIndex());
      return filter == 0 || (filter == 1 && source % 2 == 0);
    };
    const auto function = [&](const Point& point) -> Range {
      const auto& entity = point.getPolytope();
      const Index source = mesh.getGlobalIndex(entity.getDimension(), entity.getIndex());
      const auto value = [&](size_t component) -> Scalar {
        const Real real = static_cast<Real>(16 * (source + 1) + component);
        if constexpr (std::is_same_v<Scalar, Complex>)
          return {real, -real - 3};
        else
          return real;
      };
      if constexpr (std::is_same_v<Range, Scalar>)
        return value(0);
      else
      {
        Range result(fes.getVectorDimension());
        for (size_t component = 0; component < fes.getVectorDimension(); ++component)
          result(component) = value(component);
        return result;
      }
    };

    using Record = std::pair<Index, std::pair<Index, Scalar>>;
    std::vector<Record> records;
    if (region == Geometry::Region::Cells || dimension > 0)
      for (auto entity = mesh.getPolytope(sourceDimension); entity; ++entity)
      {
        const Index local = entity->getIndex();
        if (!shard.isOwned(sourceDimension, local) || !predicate(*entity))
          continue;
        if (region == Geometry::Region::Boundary && !shard.isBoundary(local))
          continue;
        if (region == Geometry::Region::Interface && !shard.isInterface(local))
          continue;
        const Index source = mesh.getGlobalIndex(sourceDimension, local);
        const auto dofs = fes.getDOFs(sourceDimension, local);
        const auto& element = fes.getFiniteElement(sourceDimension, local);
        const auto pullback = fes.getPullback({sourceDimension, local}, function);
        for (Index ordinal = 0; ordinal < static_cast<Index>(dofs.size()); ++ordinal)
        {
          records.push_back(
            {dofs[ordinal], {source, element.getLinearForm(ordinal)(pullback)}});
        }
      }
    std::vector<std::vector<Record>> gathered;
    boost::mpi::all_gather(comm, records, gathered);
    IndexMap<std::pair<Index, Scalar>> expected;
    for (const auto& rank : gathered)
    {
      for (const auto& [dof, candidate] : rank)
      {
        const auto found = expected.find(dof);
        if (found == expected.end() || candidate.first < found->second.first)
          expected[dof] = candidate;
      }
    }

    Interpolation interpolation(fes, region, predicate);
    IndexMap<Scalar> values;
    // Existing output must be replaced, including when the region is empty.
    values.emplace(std::numeric_limits<Index>::max(), Scalar(-97));
    interpolation.assemble(values, function);
    Index begin, end;
    fes.getOwnershipRange(begin, end);
    size_t expectedCount = 0;
    for (const auto& [dof, candidate] : expected)
    {
      if (begin <= dof && dof < end)
      {
        ++expectedCount;
        const auto found = values.find(dof);
        EXPECT_NE(found, values.end());
        if (found != values.end())
          EXPECT_EQ(found->second, candidate.second);
      }
    }
    EXPECT_EQ(values.size(), expectedCount);
    for (const auto& [dof, value] : values)
    {
      EXPECT_GE(dof, begin);
      EXPECT_LT(dof, end);
    }
    for (const auto& [dof, indices] : interpolation.getDOFs())
    {
      const auto [entity, ordinal] = indices;
      const auto found = expected.find(dof);
      EXPECT_NE(found, expected.end());
      if (found != expected.end())
        EXPECT_EQ(mesh.getGlobalIndex(sourceDimension, entity), found->second.first);
      EXPECT_EQ(fes.getGlobalIndex({sourceDimension, entity}, ordinal), dof);
    }
    IndexMap<Scalar> repeated;
    interpolation.assemble(repeated, function);
    EXPECT_EQ(values, repeated);
  }

  template <class FES>
  void checkRegions(const FES& fes, bool cellsOnly)
  {
    for (auto region : {Geometry::Region::Cells, Geometry::Region::Faces,
           Geometry::Region::Boundary, Geometry::Region::Interface})
    {
      if (cellsOnly && region != Geometry::Region::Cells)
        continue;
      for (int filter = 0; filter < 3; ++filter)
        check(fes, region, filter);
    }
  }

  template <class Scalar>
  void checkSpaces(const Mesh<Context::MPI>& mesh, bool allRegions)
  {
    using Vector = Math::SpatialVector<Scalar>;
    using MPIMesh = Mesh<Context::MPI>;
    checkRegions(P0<Scalar, MPIMesh>(mesh), true);
    checkRegions(P0<Vector, MPIMesh>(mesh, 2), true);
    checkRegions(P0g<Scalar, MPIMesh>(mesh), !allRegions);
    checkRegions(P0g<Vector, MPIMesh>(mesh, 2), !allRegions);
    checkRegions(P1<Scalar, MPIMesh>(mesh), !allRegions);
    checkRegions(P1<Vector, MPIMesh>(mesh, 2), !allRegions);
    checkRegions(
      H1<1, Scalar, MPIMesh>(std::integral_constant<size_t, 1>{}, mesh), !allRegions);
    checkRegions(
      H1<1, Vector, MPIMesh>(std::integral_constant<size_t, 1>{}, mesh, 2), !allRegions);
    checkRegions(
      H1<2, Scalar, MPIMesh>(std::integral_constant<size_t, 2>{}, mesh), !allRegions);
    checkRegions(
      H1<2, Vector, MPIMesh>(std::integral_constant<size_t, 2>{}, mesh, 2), !allRegions);
    checkRegions(
      H1<3, Scalar, MPIMesh>(std::integral_constant<size_t, 3>{}, mesh), !allRegions);
    checkRegions(
      H1<3, Vector, MPIMesh>(std::integral_constant<size_t, 3>{}, mesh, 2), !allRegions);
  }
}

namespace Rodin::Tests::Unit
{
  /** Empty shards remain distinct from actual point cells in zero dimension. */
  TEST(Interpolation, EmptyAndPointMeshFunctionals)
  {
    Context::MPI context(*g_environment, *g_world);
    for (const bool includePoint : {false, true})
    {
      SCOPED_TRACE(includePoint);
      Shard::Builder builder;
      builder.initialize(0, 3);
      if (includePoint && g_world->rank() == g_world->size() - 1)
        builder.vertex(0, Math::SpatialPoint{0, 0, 0}, Shard::State::Owned);
      auto mesh =
        Mesh<Context::MPI>::Builder(context).initialize(builder.finalize()).finalize();
      checkSpaces<Real>(mesh, false);
      checkSpaces<Complex>(mesh, false);
      checkRegions(P1<Real, decltype(mesh)>(mesh), false);
      checkRegions(P0g<Math::SpatialVector<Complex>, decltype(mesh)>(mesh, 2), false);
    }
  }

  /** Real/complex scalar/vector functionals and source IDs across every geometry. */
  TEST(Interpolation, OwnedFunctionalsAcrossSpacesAndGeometries)
  {
    Context::MPI context(*g_environment, *g_world);
    for (const auto geometry : geometries)
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = distribute(context, geometry);
      if (g_world->size() > 1 && g_world->rank() == 0)
        EXPECT_EQ(mesh.getShard().getVertexCount(), 0u);
      checkSpaces<Real>(mesh, false);
      checkSpaces<Complex>(mesh, false);
    }
  }

  /** Face selection includes overlap sources even when the face owner differs. */
  TEST(Interpolation, FilteredFaceBoundaryAndInterfaceFunctionals)
  {
    Context::MPI context(*g_environment, *g_world);
    for (const auto geometry : geometries)
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = distribute(context, geometry);
      P1<Real, decltype(mesh)> p1(mesh);
      P0g<Math::SpatialVector<Complex>, decltype(mesh)> constant(mesh, 2);
      H1<3, Math::SpatialVector<Complex>, decltype(mesh)> h1(
        std::integral_constant<size_t, 3>{}, mesh, 2);
      checkRegions(p1, false);
      checkRegions(constant, false);
      checkRegions(h1, false);
    }
  }

  /** Extracted skins and strict cell selections retain complete source stars. */
  TEST(Interpolation, SubmeshFunctionalsUseSubmeshLogicalIndices)
  {
    Context::MPI context(*g_environment, *g_world);
    for (const auto geometry : geometries)
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = distribute(context, geometry);
      for (const bool skin : {false, true})
      {
        SCOPED_TRACE(skin);
        const size_t dimension = mesh.getDimension() - (skin ? 1 : 0);
        SubMesh<Context::MPI>::Builder builder;
        builder.initialize(mesh);
        for (auto entity = mesh.getPolytope(dimension); entity; ++entity)
        {
          const Index local = entity->getIndex();
          if (!mesh.getShard().isOwned(dimension, local))
            continue;
          if (skin ? mesh.getShard().isBoundary(local)
                   : mesh.getGlobalIndex(dimension, local) % 2 == 0)
            builder.include(dimension, local);
        }
        auto submesh = builder.finalize();
        checkSpaces<Real>(submesh, false);
        checkSpaces<Complex>(submesh, false);
      }
    }
  }

  /** P0 has cell functionals only; unsupported face interpolation is rejected. */
  TEST(Interpolation, P0RejectsNonCellRegions)
  {
    Context::MPI context(*g_environment, *g_world);
    auto mesh = distribute(context, Polytope::Type::Triangle);
    P0<Real, decltype(mesh)> scalar(mesh);
    P0<Math::SpatialVector<Complex>, decltype(mesh)> vector(mesh, 2);
    const auto selected = [](const Polytope&) { return true; };
    for (const auto region :
      {Geometry::Region::Faces, Geometry::Region::Boundary, Geometry::Region::Interface})
    {
      EXPECT_THROW(Interpolation(scalar, region, selected), Alert::Exception);
      EXPECT_THROW(Interpolation(vector, region, selected), Alert::Exception);
    }
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator world;
  g_environment = &environment;
  g_world = &world;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
