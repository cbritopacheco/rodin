/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/GreedyPartitioner.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  void checkGreedyPartition(LocalMesh& mesh, size_t count, size_t dimension)
  {
    mesh.getConnectivity().compute(dimension, dimension);
    GreedyPartitioner partitioner(mesh);
    EXPECT_EQ(&partitioner.getMesh(), &mesh);
    partitioner.partition(count, dimension);
    ASSERT_EQ(partitioner.getCount(), count);
    const size_t n = mesh.getConnectivity().getCount(dimension);
    std::vector<size_t> sizes(count, 0);
    std::vector<size_t> assignments;
    for (size_t i = 0; i < n; ++i)
    {
      const size_t partition = partitioner.getPartition(i);
      ASSERT_LT(partition, count) << "entity=" << i;
      ++sizes[partition];
      assignments.push_back(partition);
    }
    for (size_t size : sizes)
    {
      EXPECT_LE(size, (n + count - 1) / count);
      EXPECT_GE(size, n / count);
    }
    partitioner.partition(count, dimension);
    for (size_t i = 0; i < n; ++i)
      EXPECT_EQ(partitioner.getPartition(i), assignments[i]);
  }

  /// @brief Connected cells all receive valid balanced assignments, including queued neighbors.
  TEST(GreedyPartitioner, ConnectedCells)
  {
    for (const auto type : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
    {
      auto mesh = LocalMesh::UniformGrid(type, {4, 4});
      for (size_t count : {1, 2, 3, 5})
        checkGreedyPartition(mesh, count, 2);
    }
  }

  /// @brief Disconnected components must not create more partitions than requested.
  TEST(GreedyPartitioner, DisconnectedCells)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(1)
                  .nodes(8)
                  .vertex({0})
                  .vertex({1})
                  .vertex({2})
                  .vertex({3})
                  .vertex({4})
                  .vertex({5})
                  .vertex({6})
                  .vertex({7})
                  .polytope(Polytope::Type::Segment, {{0, 1}})
                  .polytope(Polytope::Type::Segment, {{2, 3}})
                  .polytope(Polytope::Type::Segment, {{4, 5}})
                  .polytope(Polytope::Type::Segment, {{6, 7}})
                  .finalize();
    checkGreedyPartition(mesh, 2, 1);
  }

  /// @brief Excess partition requests leave empty partitions without dropping entities.
  TEST(GreedyPartitioner, MorePartitionsThanCells)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    checkGreedyPartition(mesh, 5, 1);
  }
  /// @brief Partitioning faces uses their incidence graph rather than cell indices.
  TEST(GreedyPartitioner, LowerDimensionalEntities)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    checkGreedyPartition(mesh, 3, 1);
  }

  /// @brief Empty meshes produce an empty assignment and retain the requested count.
  TEST(GreedyPartitioner, EmptyMesh)
  {
    auto mesh = LocalMesh::Builder().initialize(1).nodes(0).finalize();
    GreedyPartitioner partitioner(mesh);
    EXPECT_EQ(partitioner.getCount(), 0);
    partitioner.partition(3);
    EXPECT_EQ(partitioner.getCount(), 3);
  }

  /// @brief A zero partition request is rejected without changing an existing assignment.
  TEST(GreedyPartitioner, RejectsZeroPartitions)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    mesh.getConnectivity().compute(1, 1);
    GreedyPartitioner partitioner(mesh);
    partitioner.partition(2);
    const size_t first = partitioner.getPartition(0);
    EXPECT_THROW(partitioner.partition(0), Alert::Exception);
    EXPECT_EQ(partitioner.getCount(), 2);
    EXPECT_EQ(partitioner.getPartition(0), first);
  }

}
