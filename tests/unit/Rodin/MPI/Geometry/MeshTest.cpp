/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <array>
#include <vector>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/IO.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/QF/Centroid.h>
#include <Rodin/Variational/H1/H1Element.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;

  constexpr Polytope::Type geometries[] = {Polytope::Type::Segment,
    Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
    Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge};

  /// Distributes the cells over all ranks except root when multiple ranks exist.
  class EmptyRootPartitioner final : public Partitioner
  {
    public:
      explicit EmptyRootPartitioner(const Mesh<Context::Local>& mesh)
        : m_partitioner(mesh),
          m_count(0)
      {}

      const Mesh<Context::Local>& getMesh() const override
      {
        return m_partitioner.getMesh();
      }

      void partition(size_t count, size_t d) override
      {
        m_count = count;
        m_partitioner.partition(count > 1 ? count - 1 : 1, d);
      }

      size_t getPartition(Index index) const override
      {
        return m_partitioner.getPartition(index) + (m_count > 1 ? 1 : 0);
      }

      size_t getCount() const override
      {
        return m_count;
      }

    private:
      BalancedCompactPartitioner m_partitioner;
      size_t m_count;
  };

  Mesh<Context::MPI> distribute(Polytope::Type geometry, bool completeIncidences = true)
  {
    const Context::MPI context(*environment, *world);
    Sharder<Context::MPI> sharder(context);
    if (world->rank() == 0)
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      const Array<size_t> shape = Array<size_t>::Constant(dimension, 3);
      auto local = Mesh<Context::Local>::UniformGrid(geometry, shape);
      local.getConnectivity().compute(dimension, 0);
      local.getConnectivity().compute(dimension, dimension);
      if (completeIncidences)
      {
        for (size_t d = 0; d <= dimension; ++d)
        {
          for (size_t dp = 0; dp <= dimension; ++dp)
            local.getConnectivity().compute(d, dp);
        }
      }
      EmptyRootPartitioner partitioner(local);
      partitioner.partition(world->size(), dimension);
      sharder.shard(partitioner);
    }
    sharder.scatter(0);
    return sharder.gather(0);
  }
}

/// @brief Attribute-restricted measures recover half-unit-box volumes and boundaries.
TEST(MPI_Geometry_Mesh, AttributeRestrictedUnitBoxMeasuresWithEmptyRoot)
{
  for (const auto geometry :
    {Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const Array<size_t> shape = Array<size_t>::Constant(dimension, 3);
    const Context::MPI context(*environment, *world);
    Sharder<Context::MPI> sharder(context);
    if (world->rank() == 0)
    {
      auto local = LocalMesh::UniformGrid(geometry, shape);
      local.scale(0.5);
      local.getConnectivity().compute(dimension - 1, dimension);
      for (auto cell = local.getCell(); cell; ++cell)
      {
        const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        local.setAttribute({dimension, cell->getIndex()}, center.x() < 0.5 ? 17 : 23);
      }
      for (auto face = local.getBoundary(); face; ++face)
      {
        const Point center(*face, Polytope::Traits(face->getGeometry()).getCentroid());
        local.setAttribute({dimension - 1, face->getIndex()}, center.x() < 0.5 ? 31 : 37);
      }
      local.getConnectivity().compute(dimension, 0);
      local.getConnectivity().compute(dimension, dimension);
      EmptyRootPartitioner partitioner(local);
      partitioner.partition(world->size(), dimension);
      sharder.shard(partitioner);
    }
    sharder.scatter(0);
    auto mesh = sharder.gather(0);
    if (world->size() > 1 && world->rank() == 0)
    {
      EXPECT_EQ(mesh.getShard().getCellCount(), 0u);
    }

    for (const Attribute attribute : {17, 23})
    {
      EXPECT_NEAR(mesh.getMeasure(dimension, attribute), 0.5, 1e-13);
      EXPECT_NEAR(mesh.getMeasure(dimension, FlatSet<Attribute>{attribute}), 0.5, 1e-13);
      if (dimension == 2)
      {
        EXPECT_NEAR(mesh.getArea(attribute), 0.5, 1e-13);
        EXPECT_NEAR(mesh.getArea(FlatSet<Attribute>{attribute}), 0.5, 1e-13);
      }
      else
      {
        EXPECT_NEAR(mesh.getVolume(attribute), 0.5, 1e-13);
        EXPECT_NEAR(mesh.getVolume(FlatSet<Attribute>{attribute}), 0.5, 1e-13);
      }
    }
    EXPECT_NEAR(mesh.getMeasure(dimension, FlatSet<Attribute>{17, 23, 999}), 1, 1e-13);
    EXPECT_EQ(mesh.getMeasure(dimension, Attribute(999)), 0);
    EXPECT_EQ(mesh.getMeasure(dimension, FlatSet<Attribute>{}), 0);
    EXPECT_EQ(mesh.getMeasure(dimension, FlatSet<Attribute>{999}), 0);
    // A half of the unit square/cube's exterior has measure 2/3.
    for (const Attribute attribute : {31, 37})
    {
      EXPECT_NEAR(mesh.getPerimeter(attribute), dimension, 1e-13);
      EXPECT_NEAR(mesh.getPerimeter(FlatSet<Attribute>{attribute}), dimension, 1e-13);
    }
    EXPECT_NEAR(mesh.getPerimeter(FlatSet<Attribute>{31, 37, 999}), 2 * dimension, 1e-13);
    EXPECT_EQ(mesh.getPerimeter(Attribute(999)), 0);
    EXPECT_EQ(mesh.getPerimeter(FlatSet<Attribute>{}), 0);
  }
}

/**
 * @brief The distributed dimension is the maximum actual shard dimension.
 *
 * All seven cell geometries are partitioned with an empty root. Rank-local
 * dimension remains zero there, while distributed queries use the cell
 * dimension. Rank-only queries and copies prove that neither is collective.
 */
TEST(MPI_Geometry_Mesh, DimensionWithEmptyRootAcrossGeometries)
{
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = distribute(geometry);
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    EXPECT_EQ(mesh.getDimension(), dimension);
    if (world->size() > 1 && world->rank() == 0)
    {
      EXPECT_EQ(mesh.getShard().getDimension(), 0u);
      EXPECT_EQ(mesh.getShard().getVertexCount(), 0u);
      EXPECT_FALSE(mesh.getCell());
      EXPECT_FALSE(mesh.getFace());
      EXPECT_FALSE(mesh.getBoundary());
      EXPECT_FALSE(mesh.getInterface());
    }
    if (world->rank() == 0)
    {
      const Context::MPI context(*environment, *world);
      Mesh copied(mesh);
      Mesh moved(std::move(copied));
      Mesh<Context::MPI> assigned(context);
      assigned = mesh;
      Mesh<Context::MPI> moveAssigned(context);
      moveAssigned = std::move(assigned);
      EXPECT_EQ(moved.getDimension(), dimension);
      EXPECT_EQ(moveAssigned.getDimension(), dimension);
      moveAssigned.flush();
      EXPECT_EQ(moveAssigned.getDimension(), dimension);
    }
    world->barrier();
  }
}

/// @brief Invalid grid shapes are rejected consistently on every rank.
TEST(MPI_Geometry_Mesh, UniformGridRejectsInvalidShapes)
{
  const Context::MPI context(*environment, *world);
  EXPECT_THROW(Mesh<Context::MPI>::UniformGrid(
                 context, Polytope::Type::Point, Array<size_t>::Constant(1, 2)),
    Alert::Exception);
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    EXPECT_THROW(Mesh<Context::MPI>::UniformGrid(
                   context, geometry, Array<size_t>::Constant(dimension + 1, 2)),
      Alert::Exception);
    for (size_t component = 0; component < dimension; ++component)
    {
      for (const size_t count : {0u, 1u})
      {
        Array<size_t> shape = Array<size_t>::Constant(dimension, 2);
        shape[component] = count;
        EXPECT_THROW(
          Mesh<Context::MPI>::UniformGrid(context, geometry, shape), Alert::Exception);
      }
    }
  }
}

/// @brief A uniform point grid has one global owner and valid ghost replicas.
TEST(MPI_Geometry_Mesh, UniformPointGridOwnership)
{
  const Context::MPI context(*environment, *world);
  const auto mesh =
    Mesh<Context::MPI>::UniformGrid(context, Polytope::Type::Point, Array<size_t>(0));
  EXPECT_EQ(mesh.getDimension(), 0u);
  EXPECT_EQ(mesh.getSpaceDimension(), 0u);
  EXPECT_EQ(mesh.getVertexCount(), 1u);
  const auto& shard = mesh.getShard();
  ASSERT_EQ(shard.getVertexCount(), 1u);
  EXPECT_EQ(shard.isOwned(0, 0), world->rank() == 0);
  EXPECT_EQ(shard.getVertexCoordinates(0).size(), 0);
}

/// @brief Uniform-grid construction reports the same dimension on empty ranks.
TEST(MPI_Geometry_Mesh, UniformGridDimensionAcrossGeometries)
{
  const Context::MPI context(*environment, *world);
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = Mesh<Context::MPI>::UniformGrid(
      context, geometry, Array<size_t>::Constant(dimension, 2));
    EXPECT_EQ(mesh.getDimension(), dimension);
  }
}

/// @brief Empty meshes and meshes containing a single point both have dimension zero.
TEST(MPI_Geometry_Mesh, EmptyAndPointDimensions)
{
  const Context::MPI context(*environment, *world);
  Mesh<Context::MPI> defaultEmpty(context);
  EXPECT_EQ(defaultEmpty.getDimension(), 0u);
  for (const bool includePoint : {false, true})
  {
    Shard::Builder builder;
    builder.initialize(0, 3);
    if (includePoint && world->rank() == world->size() - 1)
      builder.vertex(0, Math::SpatialPoint{0, 0, 0}, Shard::State::Owned);
    auto mesh =
      Mesh<Context::MPI>::Builder(context).initialize(builder.finalize()).finalize();
    EXPECT_EQ(mesh.getDimension(), 0u);
    EXPECT_EQ(mesh.getShard().getDimension(), 0u);
    EXPECT_EQ(mesh.getVertexCount(), includePoint ? 1u : 0u);
    EXPECT_FALSE(mesh.getFace());
    EXPECT_FALSE(mesh.getBoundary());
  }
}

/// @brief Submesh extraction establishes its own dimension, including an empty selection.
TEST(MPI_Geometry_Mesh, SubMeshDimensionsWithEmptyRootAcrossGeometries)
{
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = distribute(geometry);
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (auto face = mesh.getBoundary(); face; ++face)
      builder.include(dimension - 1, face->getIndex());
    auto boundary = builder.finalize();
    EXPECT_EQ(boundary.getDimension(), dimension - 1);
    if (world->size() > 1 && world->rank() == 0)
    {
      EXPECT_EQ(boundary.getShard().getDimension(), 0u);
    }
    if (world->rank() == 0)
    {
      Mesh<Context::MPI> baseCopy(boundary);
      EXPECT_EQ(baseCopy.getDimension(), dimension - 1);
      SubMesh copy(boundary);
      SubMesh moved(std::move(copy));
      EXPECT_EQ(moved.getDimension(), dimension - 1);
    }
    SubMesh<Context::MPI>::Builder emptyBuilder;
    auto empty = emptyBuilder.initialize(mesh).finalize();
    EXPECT_EQ(empty.getDimension(), 0u);
    EXPECT_EQ(empty.getVertexCount(), 0u);
  }
}

/**
 * @brief Empty ranks participate in reconciliation and contribute zero owned entities.
 *
 * After local discovery, global entity counts must equal those of an independent
 * local grid. Repeating reconciliation preserves the distributed index map.
 */
TEST(MPI_Geometry_Mesh, ReconciliationWithEmptyRootAcrossGeometries)
{
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = distribute(geometry, false);
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto reference =
      Mesh<Context::Local>::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
    for (size_t d = 1; d < dimension; ++d)
    {
      SCOPED_TRACE(d);
      reference.getConnectivity().compute(d, dimension);
      mesh.getConnectivity().compute(d, dimension);
      mesh.reconcile(d);
      EXPECT_EQ(mesh.getPolytopeCount(d), reference.getPolytopeCount(d));
      const auto indices = mesh.getShard().getPolytopeMap(d).left;
      mesh.reconcile(d);
      EXPECT_EQ(mesh.getShard().getPolytopeMap(d).left, indices);
      EXPECT_EQ(mesh.getDimension(), dimension);
    }
  }
}

/// @brief Non-strict finite caps preserve fully propagated IDs and ownership metadata.
TEST(MPI_Geometry_Mesh, NonStrictRoundCapsPreserveConvergedEntities)
{
  boost::mpi::communicator self(MPI_COMM_SELF, boost::mpi::comm_attach);
  for (const bool independent : {true, false})
  {
    SCOPED_TRACE(independent);
    const auto& communicator = independent ? self : *world;
    const Context::MPI context(*environment, communicator);
    const size_t cap = independent ? 0 : static_cast<size_t>(communicator.size());
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
           Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      const Array<size_t> shape = Array<size_t>::Constant(dimension, 4);
      auto mesh = Mesh<Context::MPI>::UniformGrid(context, geometry, shape);
      auto reference = LocalMesh::UniformGrid(geometry, shape);
      for (size_t d = 1; d < dimension; ++d)
      {
        SCOPED_TRACE(d);
        mesh.getConnectivity().compute(dimension, d);
        mesh.getConnectivity().compute(d, 0);
        mesh.getConnectivity().compute(d, dimension);
        reference.getConnectivity().compute(d, dimension);
        mesh.reconcile(d);
        const auto& shard = mesh.getShard();
        const auto indices = shard.getPolytopeMap(d).left;
        const auto owners = shard.getOwner(d);
        const auto halo = shard.getHalo(d);
        std::vector<std::array<bool, 3>> states;
        for (Index i = 0; i < shard.getPolytopeCount(d); ++i)
        {
          states.push_back(
            {shard.isOwned(d, i), shard.isShared(d, i), shard.isGhost(d, i)});
        }
        Mesh<Context::MPI>::ReconcileOptions options;
        options.maxOwnerRounds = cap;
        options.maxGidRounds = cap;
        options.globalCheckPeriod = cap + 1;
        options.strictRoundCap = false;
        // A rank graph has diameter at most size-1. The larger check period
        // makes the first epoch use the finite cap; SELF needs no exchange.
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          SCOPED_TRACE(repeat);
          mesh.reconcile(d, options);
          EXPECT_EQ(mesh.getPolytopeCount(d), reference.getPolytopeCount(d));
          EXPECT_EQ(shard.getPolytopeMap(d).left, indices);
          EXPECT_EQ(shard.getOwner(d), owners);
          EXPECT_EQ(shard.getHalo(d), halo);
          EXPECT_EQ(shard.getPolytopeCount(d), states.size());
          if (shard.getPolytopeCount(d) == states.size())
          {
            for (Index i = 0; i < states.size(); ++i)
            {
              EXPECT_EQ(shard.isOwned(d, i), states[i][0]);
              EXPECT_EQ(shard.isShared(d, i), states[i][1]);
              EXPECT_EQ(shard.isGhost(d, i), states[i][2]);
            }
          }
        }
      }
    }
  }
}

/// @brief Both HDF5 entry points restore collective dimension with an empty root shard.
TEST(MPI_Geometry_Mesh, HDF5DimensionWithEmptyRootAcrossGeometries)
{
  const Context::MPI context(*environment, *world);
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = distribute(geometry);
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const auto path = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-mpi-dimension-%%%%-%%%%-%%%%.h5");
    mesh.save(path, IO::FileFormat::HDF5);
    Mesh<Context::MPI> loaded(context);
    loaded.load(path, IO::FileFormat::HDF5);
    EXPECT_EQ(loaded.getDimension(), dimension);
    EXPECT_EQ(loaded.getShard().getDimension(), mesh.getShard().getDimension());
    Mesh<Context::MPI> direct(context);
    IO::MeshLoader<IO::FileFormat::HDF5, Context::MPI> loader(direct);
    loader.load(path);
    EXPECT_EQ(direct.getDimension(), dimension);
    for (size_t d = 1; d < dimension; ++d)
    {
      direct.getConnectivity().compute(d, dimension);
      direct.reconcile(d);
      EXPECT_EQ(direct.getPolytopeCount(d), mesh.getPolytopeCount(d));
    }
    boost::filesystem::remove(path);
  }
}

/// @brief Already reconciled vertex and cell maps retain ownership and coordinates.
TEST(MPI_Geometry_Mesh, ReconcileVerticesAndCellsPreservesStateWithEmptyRoot)
{
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = distribute(geometry);
    const auto coordinates = mesh.getShard().getVertices().getMatrix().eval();
    for (const size_t dimension : {size_t(0), mesh.getDimension()})
    {
      SCOPED_TRACE(dimension);
      const auto indices = mesh.getShard().getPolytopeMap(dimension).left;
      const auto owners = mesh.getShard().getOwner(dimension);
      const auto halo = mesh.getShard().getHalo(dimension);
      for (size_t repeat = 0; repeat < 2; ++repeat)
      {
        EXPECT_EQ(&mesh.reconcile(dimension), &mesh);
        EXPECT_EQ(mesh.getShard().getPolytopeMap(dimension).left, indices);
        EXPECT_EQ(mesh.getShard().getOwner(dimension), owners);
        EXPECT_EQ(mesh.getShard().getHalo(dimension), halo);
        EXPECT_EQ(mesh.getShard().getVertices().getMatrix(), coordinates);
      }
    }
  }
}

/// @brief Cached MPI quadratures reuse formulas and track scaled physical measures.
TEST(MPI_Geometry_Mesh, QuadratureCacheTracksScalingWithEmptyRoot)
{
  // Roundoff budget for summing mapped quadrature weights across the grid.
  constexpr Real tolerance = 1e-12;
  for (const auto geometry : geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = distribute(geometry);
    const size_t dimension = mesh.getDimension();
    const auto& first = QF::PolytopeQuadratureFormula::get(2, geometry);
    const auto& second = QF::PolytopeQuadratureFormula::get(3, geometry);
    // distribute() creates three vertices per axis, spanning [0, 2].
    Real expectedVolume = 1;
    for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
      expectedVolume *= 2;
    for (size_t phase = 0; phase < 4; ++phase)
    {
      SCOPED_TRACE(phase);
      Real localVolume = 0;
      const auto& shard = mesh.getShard();
      for (Index cell = 0; cell < shard.getPolytopeCount(dimension); ++cell)
      {
        const auto& quadrature = mesh.getQuadrature(dimension, cell, first);
        const auto& other = mesh.getQuadrature(dimension, cell, second);
        EXPECT_NE(&quadrature, &other);
        EXPECT_EQ(&mesh.getQuadrature(dimension, cell, first), &quadrature);
        EXPECT_EQ(&quadrature.getQuadratureFormula(), &first);
        EXPECT_EQ(quadrature.getSize(), first.getSize());
        for (size_t point = 0; point < quadrature.getSize(); ++point)
        {
          const auto& mapped = quadrature.getPoint(point);
          EXPECT_EQ(&mapped.getPolytope().getMesh(), &mesh);
          if (shard.isOwned(dimension, cell))
            localVolume += first.getWeight(point) * mapped.getDistortion();
        }
      }
      const Real volume = boost::mpi::all_reduce(*world, localVolume, std::plus<Real>());
      EXPECT_NEAR(volume, expectedVolume, tolerance * expectedVolume);
      if (phase < 3)
      {
        if (phase == 0)
          EXPECT_EQ(&mesh.scale(2), &mesh);
        else
        {
          for (Index vertex = 0; vertex < shard.getVertexCount(); ++vertex)
          {
            const Math::SpatialPoint coordinates = mesh.getVertexCoordinates(vertex);
            if (phase == 1)
              EXPECT_EQ(&mesh.setVertexCoordinates(vertex, Real(2) * coordinates), &mesh);
            else
            {
              for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
              {
                EXPECT_EQ(&mesh.setVertexCoordinates(
                            vertex, Real(2) * coordinates[coordinate], coordinate),
                  &mesh);
              }
            }
          }
        }
        for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
          expectedVolume *= 2;
      }
    }
  }
}

/// @brief Replacing a cell map invalidates cached physical quadrature coordinates.
TEST(MPI_Geometry_Mesh, QuadratureCacheTracksReplacementTransformationWithEmptyRoot)
{
  auto mesh = distribute(Polytope::Type::Segment);
  const Variational::RealH1Element<2> element(Polytope::Type::Segment);
  const auto& rule = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Segment);
  // The quadratic map fixes both endpoints and has positive derivative on [0, 1].
  constexpr Real curvature = 0.25;
  constexpr Real tolerance = 1e-12;
  size_t localCount = 0;
  for (Index cell = 0; cell < mesh.getShard().getCellCount(); ++cell)
  {
    const auto& original = mesh.getPolytopeTransformation(1, cell);
    Math::SpatialPoint left, right;
    original.transform(left, Math::SpatialPoint{0});
    original.transform(right, Math::SpatialPoint{1});
    const Real center =
      mesh.getQuadrature(1, cell, rule).getPoint(0).getPhysicalCoordinates()[0];
    PointCloud nodes(1, element.getCount());
    for (size_t node = 0; node < element.getCount(); ++node)
    {
      const Real t = element.getNode(node)[0];
      nodes(0, node) = left[0] + (right[0] - left[0]) * (t + curvature * t * (1 - t));
    }
    EXPECT_EQ(&mesh.setPolytopeTransformation({1, cell},
                new ParametricTransformation<Variational::RealH1Element<2>>(
                  std::move(nodes), element)),
      &mesh);
    const auto& quadrature = mesh.getQuadrature(1, cell, rule);
    EXPECT_EQ(&mesh.getQuadrature(1, cell, rule), &quadrature);
    EXPECT_EQ(&quadrature.getPoint(0).getPolytope().getMesh(), &mesh);
    EXPECT_NEAR(quadrature.getPoint(0).getPhysicalCoordinates()[0],
      center + (right[0] - left[0]) * curvature / 4, tolerance);
    ++localCount;
  }
  EXPECT_GT(boost::mpi::all_reduce(*world, localCount, std::plus<size_t>()), 0u);
}

int main(int argc, char** argv)
{
  boost::mpi::environment mpiEnvironment(argc, argv);
  boost::mpi::communicator mpiWorld;
  environment = &mpiEnvironment;
  world = &mpiWorld;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
