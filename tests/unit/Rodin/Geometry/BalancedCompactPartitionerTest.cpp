/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <queue>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /**
   * Partitions the cells of a connected mesh and checks that every cell has a
   * partition, that the sizes differ by at most one and that every partition
   * is connected.
   */
  void checkBalancedCompactPartition(LocalMesh& mesh, size_t count)
  {
    const size_t d = mesh.getDimension();
    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(count, d);
    ASSERT_EQ(partitioner.getCount(), count);
    const auto& conn = mesh.getConnectivity();
    const size_t n = conn.getCount(d);
    std::vector<size_t> sizes(count, 0);
    std::vector<size_t> assignments(n);
    for (Index i = 0; i < n; i++)
    {
      assignments[i] = partitioner.getPartition(i);
      ASSERT_LT(assignments[i], count) << "cell=" << i;
      sizes[assignments[i]]++;
    }
    for (size_t p = 0; p < count; p++)
    {
      EXPECT_LE(sizes[p], (n + count - 1) / count) << "partition=" << p;
      EXPECT_GE(sizes[p], n / count) << "partition=" << p;
    }

    const auto& adjacency = conn.getIncidence(d, d);
    std::vector<bool> reached(n, false);
    for (size_t p = 0; p < count; p++)
    {
      Index start = n;
      for (Index i = 0; i < n && start == n; i++)
      {
        if (assignments[i] == p)
          start = i;
      }
      if (start == n)
        continue;
      size_t found = 0;
      std::queue<Index> queue;
      queue.push(start);
      reached[start] = true;
      while (!queue.empty())
      {
        const Index i = queue.front();
        queue.pop();
        found++;
        for (const Index j : adjacency[i])
        {
          if (!reached[j] && assignments[j] == p)
          {
            reached[j] = true;
            queue.push(j);
          }
        }
      }
      EXPECT_EQ(found, sizes[p]) << "partition " << p << " is not connected";
    }

    partitioner.partition(count, d);
    for (Index i = 0; i < n; i++)
      EXPECT_EQ(partitioner.getPartition(i), assignments[i]);
  }

  /// @brief Triangle and quadrilateral grids split into equal connected parts.
  TEST(BalancedCompactPartitioner, Surface)
  {
    for (const auto type : { Polytope::Type::Triangle, Polytope::Type::Quadrilateral })
    {
      auto mesh = LocalMesh::UniformGrid(type, { 32, 24 });
      mesh.getConnectivity().compute(2, 2);
      for (size_t count : { 1, 2, 3, 4, 7, 8, 16 })
        checkBalancedCompactPartition(mesh, count);
    }
  }

  /// @brief Tetrahedral grids split into equal connected parts.
  TEST(BalancedCompactPartitioner, Volume)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 10, 8, 6 });
    mesh.getConnectivity().compute(3, 3);
    for (size_t count : { 2, 3, 4, 8, 12 })
      checkBalancedCompactPartition(mesh, count);
  }

  /// @brief Every cell of a mesh with seedless components gets a partition.
  TEST(BalancedCompactPartitioner, DisconnectedCells)
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
    mesh.getConnectivity().compute(1, 1);
    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2, 1);
    std::vector<size_t> sizes(2, 0);
    for (Index i = 0; i < 4; i++)
    {
      ASSERT_LT(partitioner.getPartition(i), 2u);
      sizes[partitioner.getPartition(i)]++;
    }
    EXPECT_EQ(sizes[0] + sizes[1], 4u);
    EXPECT_GE(sizes[0], 1u);
    EXPECT_GE(sizes[1], 1u);
  }

  /// @brief More partitions than cells leaves the extra partitions empty.
  TEST(BalancedCompactPartitioner, MorePartitionsThanCells)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 4 });
    mesh.getConnectivity().compute(1, 1);
    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(5, 1);
    std::vector<size_t> sizes(5, 0);
    for (Index i = 0; i < 3; i++)
    {
      ASSERT_LT(partitioner.getPartition(i), 5u);
      sizes[partitioner.getPartition(i)]++;
    }
    for (size_t p = 0; p < 5; p++)
      EXPECT_LE(sizes[p], 1u);
  }
}
