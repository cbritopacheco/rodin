/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Scotch partition coverage, lifecycle and strategy ownership checks. */
#include <gtest/gtest.h>
#include <limits>
#include "Rodin/Scotch/MeshPartitioner.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Every entity receives a valid assignment and each requested part is populated.
  void checkScotchPartition(
    Scotch::Partitioner& partitioner, size_t count, size_t dimension)
  {
    partitioner.partition(count, dimension);
    ASSERT_EQ(partitioner.getCount(), count);
    const size_t n = partitioner.getMesh().getConnectivity().getCount(dimension);
    std::vector<size_t> sizes(count, 0);
    for (Index i = 0; i < n; ++i)
    {
      ASSERT_LT(partitioner.getPartition(i), count);
      ++sizes[partitioner[i]];
    }
    // Scotch's heuristic need not match a particular numbering or exact quota.
    for (size_t size : sizes)
      EXPECT_GT(size, 0);
  }

  /// @brief Construction and destruction do not require a prior partition operation.
  TEST(ScotchPartitioner, UnusedLifetime)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    Scotch::Partitioner partitioner(mesh);
    EXPECT_EQ(&partitioner.getMesh(), &mesh);
    EXPECT_EQ(partitioner.getCount(), 0);
  }

  /// @brief Successive requests rebuild the graph and replace the assignment.
  TEST(ScotchPartitioner, RepeatedPartitions)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
    mesh.getConnectivity().compute(2, 2);
    Scotch::Partitioner partitioner(mesh);
    for (size_t count : {1, 2, 4, 2})
      checkScotchPartition(partitioner, count, 2);
    partitioner.partition(1);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
      EXPECT_EQ(partitioner.getPartition(i), 0);
  }

  /// @brief An empty mesh can be partitioned and destroyed safely.
  TEST(ScotchPartitioner, EmptyMesh)
  {
    auto mesh = LocalMesh::Builder().initialize(1).nodes(0).finalize();
    Scotch::Partitioner partitioner(mesh);
    partitioner.partition(3);
    EXPECT_EQ(partitioner.getCount(), 3);
  }

  /// @brief A zero-partition request preserves the previous successful assignment.
  TEST(ScotchPartitioner, RejectsZeroPartitions)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {9});
    mesh.getConnectivity().compute(1, 1);
    Scotch::Partitioner partitioner(mesh);
    partitioner.partition(2);
    const auto first = partitioner.getPartition(0);
    EXPECT_THROW(partitioner.partition(0), Alert::Exception);
    EXPECT_EQ(partitioner.getCount(), 2);
    EXPECT_EQ(partitioner.getPartition(0), first);
  }

  /// @brief The requested face dimension determines the graph and assignment cardinality.
  TEST(ScotchPartitioner, LowerDimension)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
    mesh.getConnectivity().compute(1, 1);
    Scotch::Partitioner partitioner(mesh);
    checkScotchPartition(partitioner, 3, 1);
  }

  /// @brief Strategy contents are copied independently of the caller's Scotch allocation.
  TEST(ScotchPartitioner, StrategyLifetime)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
    mesh.getConnectivity().compute(2, 2);
    Scotch::Partitioner first(mesh), second(mesh);
    SCOTCH_Strat strategy;
    ASSERT_EQ(SCOTCH_stratInit(&strategy), 0);
    ASSERT_EQ(SCOTCH_stratGraphMapBuild(&strategy, SCOTCH_STRATQUALITY, 2, 0.05), 0);
    EXPECT_EQ(&first.setStrategy(strategy), &first);
    second.setStrategy(strategy);
    SCOTCH_stratExit(&strategy);
    checkScotchPartition(first, 2, 2);
    checkScotchPartition(second, 2, 2);
    auto copied = first;
    checkScotchPartition(copied, 2, 2);
  }
  /// @brief Requests outside Scotch's integer range or the mesh dimensions fail safely.
  TEST(ScotchPartitioner, InvalidArguments)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {9});
    mesh.getConnectivity().compute(1, 1);
    Scotch::Partitioner partitioner(mesh);
    EXPECT_THROW(
      partitioner.partition(size_t(std::numeric_limits<SCOTCH_Num>::max()) + 1),
      Alert::Exception);
    EXPECT_THROW(partitioner.partition(2, 2), Alert::Exception);
    partitioner.partition(2);
    EXPECT_THROW(partitioner.getPartition(mesh.getCellCount()), std::out_of_range);
  }

}
