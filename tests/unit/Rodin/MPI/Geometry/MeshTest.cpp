/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/IO.h>

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
        for (size_t d = 0; d <= dimension; ++d)
          for (size_t dp = 0; dp <= dimension; ++dp)
            local.getConnectivity().compute(d, dp);
      EmptyRootPartitioner partitioner(local);
      partitioner.partition(world->size(), dimension);
      sharder.shard(partitioner);
    }
    sharder.scatter(0);
    return sharder.gather(0);
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

int main(int argc, char** argv)
{
  boost::mpi::environment mpiEnvironment(argc, argv);
  boost::mpi::communicator mpiWorld;
  environment = &mpiEnvironment;
  world = &mpiWorld;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
