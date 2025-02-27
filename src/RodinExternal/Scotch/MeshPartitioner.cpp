#include "Exception.h"
#include "MeshPartitioner.h"

#include <scotch.h>

namespace Rodin::External::Scotch
{
  const MeshPartitioner::MeshType& MeshPartitioner::getMesh() const
  {
    return m_mesh.get();
  }

  void MeshPartitioner::partition(size_t numPartitions, size_t d)
  {
    const auto& mesh = getMesh();
    const auto& conn = m_mesh.get().getConnectivity();
    const size_t n = conn.getCount(d);
    if (n == 0)
      return;
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d)
    const auto& dual = conn.getIncidence(d, d);
    std::vector<int> verttab(n + 1, 0);
    std::vector<int> edgetab;
    for (size_t i = 0; i < n; i++)
    {
      std::vector<int> neighbors;
      // Each dual[i] is a set-like container of neighboring entity indices.
      for (const auto& nb : dual[i])
      {
        if (nb != i) // Ignore self-links if any
          neighbors.push_back(static_cast<Integer>(nb));
      }
      verttab[i + 1] = verttab[i] + static_cast<int>(neighbors.size());
      edgetab.insert(edgetab.end(), neighbors.begin(), neighbors.end());
    }
    m_partition.resize(n, -1);
    SCOTCH_Graph graph;
    if (SCOTCH_graphInit(&graph) != 0)
      Exception("SCOTCH_graphInit failed", *this, __func__).raise();
    const int baseval = 0;  // using 0-based indexing
    if (SCOTCH_graphBuild(&graph,
                          baseval,
                          n,
                          verttab.data(),
                          nullptr,      // vendtab (optional, not needed with baseval=0)
                          nullptr,      // velotab (vertex weights, not used)
                          0,            // no vertex weights provided
                          edgetab.size(),
                          edgetab.data(),
                          nullptr       // edlotab (edge weights, not used)
                          ) != 0)
    {
      SCOTCH_graphExit(&graph);
      Exception("SCOTCH_graphBuild failed", *this, __func__).raise();
    }
    // Compute the partitioning (mapping) into numPartitions parts.
    if (SCOTCH_graphMap(&graph, static_cast<int>(numPartitions), m_partition.data()) != 0)
    {
      SCOTCH_graphExit(&graph);
      throw std::runtime_error("SCOTCH_graphMap failed");
    }
    SCOTCH_graphExit(&graph);
    // Optionally, print partitioning result for debug purposes.
    std::cout << "Partitioning result (entity dimension " << d << "):\n";
    for (size_t i = 0; i < n; i++)
    {
      std::cout << "  Entity " << i << " -> partition " << m_partition[i] << "\n";
    }
  }

  size_t MeshPartitioner::getPartition(Index index) const
  {
    return m_partition[index];
  }
}
