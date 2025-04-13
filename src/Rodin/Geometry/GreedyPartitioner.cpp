#include <queue>

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
    const auto& mesh = getMesh();
    const auto& conn = mesh.getConnectivity();
    const size_t n = conn.getCount(d);
    if (n == 0)
      return;
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d)
    const size_t maxPartitionSize = n / numPartitions;
    size_t partitionId = 0;
    std::vector<bool> visited(n, false);
    m_partition.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
      if (visited[i])
        continue;
      std::queue<size_t> q;
      std::vector<size_t> currentGroup;
      q.push(i);
      visited[i] = true;
      while (!q.empty() && currentGroup.size() < maxPartitionSize)
      {
        size_t current = q.front();
        q.pop();
        currentGroup.push_back(current);
        auto currentPolytope = mesh.getPolytope(mesh.getDimension(), current);
        for (auto adj = currentPolytope->getAdjacent(); !adj.end(); ++adj)
        {
          size_t adjIndex = adj->getIndex();
          if (!visited[adjIndex])
          {
            visited[adjIndex] = true;
            q.push(adjIndex);
          }
        }
      }
      for (size_t idx : currentGroup)
        m_partition[idx] = partitionId;
      partitionId++;
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
