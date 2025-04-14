#include "Exception.h"
#include "MeshPartitioner.h"

namespace Rodin::External::Scotch
{
  Partitioner::Partitioner(const MeshType& mesh)
    : m_mesh(mesh)
  {}

  Partitioner::~Partitioner()
  {
    SCOTCH_stratExit(&m_strat.value());
    SCOTCH_graphExit(&m_graph);
  }

  Partitioner& Partitioner::setStrategy(const SCOTCH_Strat& strat)
  {
    m_strat = strat;
    return *this;
  }

  const Partitioner::MeshType& Partitioner::getMesh() const
  {
    return m_mesh.get();
  }

  void Partitioner::partition(size_t numPartitions, size_t d)
  {
    m_numPartitions = numPartitions;
    const auto& mesh = getMesh();
    const auto& conn = m_mesh.get().getConnectivity();
    const size_t n = conn.getCount(d);
    if (n == 0)
      return;
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d)
    const auto& dual = conn.getIncidence(d, d);
    std::vector<SCOTCH_Num> verttab(n + 1, 0);
    std::vector<SCOTCH_Num> edgetab;
    for (size_t i = 0; i < n; i++)
    {
      std::vector<SCOTCH_Num> neighbors;
      // Each dual[i] is a set-like container of neighboring entity indices.
      for (const auto& nb : dual[i])
      {
        if (nb != i) // Ignore self-links if any
          neighbors.push_back(static_cast<SCOTCH_Num>(nb));
      }
      verttab[i + 1] = verttab[i] + static_cast<SCOTCH_Num>(neighbors.size());
      edgetab.insert(edgetab.end(), neighbors.begin(), neighbors.end());
    }
    m_partition.resize(n, -1);
    if (SCOTCH_graphInit(&m_graph) != 0)
      Exception("SCOTCH_graphInit failed", *this, __func__).raise();
    constexpr SCOTCH_Num baseval = 0;  // using 0-based indexing
    if (SCOTCH_graphBuild(&m_graph,
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
      SCOTCH_graphExit(&m_graph);
      Exception("SCOTCH_graphBuild failed", *this, __func__).raise();
    }

    if (!m_strat.has_value())
    {
      m_strat.emplace();
      if (SCOTCH_stratInit(&m_strat.value()) != 0)
      {
        SCOTCH_graphExit(&m_graph);
        Exception("SCOTCH_stratInit failed", *this, __func__).raise();
      }
    }
    auto& strat = m_strat.value();

    // Compute the partitioning (mapping) into 'numPartitions' parts.
    if (SCOTCH_graphPart(&m_graph, numPartitions, &strat, m_partition.data()) != 0)
    {
      SCOTCH_stratExit(&strat);
      SCOTCH_graphExit(&m_graph);
      Exception("SCOTCH_graphMap failed", *this, __func__).raise();
    }
  }

  size_t Partitioner::getPartition(Index index) const
  {
    return m_partition[index];
  }
}
