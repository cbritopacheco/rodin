/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <queue>

#include "Rodin/Alert/MemberFunctionException.h"
#include "GreedyPartitioner.h"

namespace Rodin::Geometry
{
  GreedyPartitioner::GreedyPartitioner(const MeshType& mesh)
    : m_mesh(mesh)
  {}

  const Mesh<Context::Local>& GreedyPartitioner::getMesh() const
  {
    return m_mesh.get();
  }

  void GreedyPartitioner::partition(size_t numPartitions, size_t d)
  {
    if (numPartitions == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "The number of partitions must be positive." << Alert::Raise;
    }
    m_count = numPartitions;
    const auto& mesh = getMesh();
    const auto& conn = mesh.getConnectivity();
    const size_t n = conn.getCount(d);
    m_partition.resize(n);
    if (n == 0)
      return;
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d)
    const auto& adjacency = conn.getIncidence(d, d);
    const size_t baseSize = n / numPartitions;
    const size_t remainder = n % numPartitions;
    size_t partitionId = 0;
    size_t partitionSize = 0;
    std::vector<bool> visited(n, false);
    std::queue<size_t> queue;
    for (size_t seed = 0; seed < n; ++seed)
    {
      if (visited[seed])
        continue;
      queue.push(seed);
      visited[seed] = true;
      while (!queue.empty())
      {
        const size_t current = queue.front();
        queue.pop();
        m_partition[current] = partitionId;
        ++partitionSize;
        if (partitionSize == baseSize + (partitionId < remainder))
        {
          ++partitionId;
          partitionSize = 0;
        }
        for (const Index neighbor : adjacency[current])
        {
          if (!visited[neighbor])
          {
            visited[neighbor] = true;
            queue.push(neighbor);
          }
        }
      }
    }
  }

  size_t GreedyPartitioner::getPartition(Index index) const
  {
    return m_partition[index];
  }

  size_t GreedyPartitioner::getCount() const
  {
    return m_count;
  }
}
