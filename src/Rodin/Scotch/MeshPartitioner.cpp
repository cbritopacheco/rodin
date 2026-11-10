/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdio>
#include <limits>
#include <memory>

#include "Exception.h"
#include "MeshPartitioner.h"

namespace Rodin::Scotch
{
  Partitioner::Partitioner(const MeshType& mesh)
    : m_numPartitions(0),
      m_mesh(mesh)
  {}

  Partitioner::~Partitioner() = default;

  Partitioner& Partitioner::setStrategy(const SCOTCH_Strat& strat)
  {
    // Scotch strategies own an opaque allocation: copying the C struct aliases it.
    // Serialize the mapping strategy so this partitioner owns an independent value.
    std::unique_ptr<FILE, decltype(&std::fclose)> stream(std::tmpfile(), &std::fclose);
    if (!stream || SCOTCH_stratSave(&strat, stream.get()) != 0 ||
      std::fseek(stream.get(), 0, SEEK_SET) != 0)
      Exception("SCOTCH strategy serialization failed", *this, __func__).raise();
    std::string strategy;
    constexpr size_t BufferSize = 4096;
    char buffer[BufferSize];
    while (const size_t count = std::fread(buffer, 1, BufferSize, stream.get()))
      strategy.append(buffer, count);
    if (std::ferror(stream.get()))
      Exception("SCOTCH strategy read failed", *this, __func__).raise();
    m_strategy = std::move(strategy);
    return *this;
  }

  const Partitioner::MeshType& Partitioner::getMesh() const
  {
    return m_mesh.get();
  }

  void Partitioner::partition(size_t numPartitions, size_t d)
  {
    constexpr size_t MaxScotchCount = std::numeric_limits<SCOTCH_Num>::max();
    if (numPartitions == 0 || numPartitions > MaxScotchCount)
    {
      Exception(
        "Partition count is outside Scotch's positive integer range", *this, __func__)
        .raise();
    }
    const auto& mesh = getMesh();
    if (d > mesh.getDimension())
    {
      Exception("Partition dimension exceeds the mesh dimension", *this, __func__)
        .raise();
    }
    const auto& conn = mesh.getConnectivity();
    const size_t n = conn.getCount(d);
    if (n == 0)
    {
      m_partition.clear();
      m_numPartitions = numPartitions;
      return;
    }
    if (n > MaxScotchCount)
      Exception("Entity count exceeds Scotch's integer range", *this, __func__).raise();
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d);
    const auto& dual = conn.getIncidence(d, d);
    std::vector<SCOTCH_Num> verttab(n + 1, 0);
    std::vector<SCOTCH_Num> edgetab;
    for (size_t i = 0; i < n; ++i)
    {
      for (const auto neighbor : dual[i])
      {
        if (neighbor != i)
        {
          if (edgetab.size() == MaxScotchCount)
          {
            Exception("Edge count exceeds Scotch's integer range", *this, __func__)
              .raise();
          }
          edgetab.push_back(static_cast<SCOTCH_Num>(neighbor));
        }
      }
      verttab[i + 1] = static_cast<SCOTCH_Num>(edgetab.size());
    }
    SCOTCH_Graph graph;
    if (SCOTCH_graphInit(&graph) != 0)
      Exception("SCOTCH_graphInit failed", *this, __func__).raise();
    std::unique_ptr<SCOTCH_Graph, decltype(&SCOTCH_graphExit)> graphCleanup(
      &graph, &SCOTCH_graphExit);
    if (SCOTCH_graphBuild(&graph, 0, static_cast<SCOTCH_Num>(n), verttab.data(), nullptr,
          nullptr, nullptr, static_cast<SCOTCH_Num>(edgetab.size()), edgetab.data(),
          nullptr) != 0)
      Exception("SCOTCH_graphBuild failed", *this, __func__).raise();

    SCOTCH_Strat strategy;
    if (SCOTCH_stratInit(&strategy) != 0)
      Exception("SCOTCH_stratInit failed", *this, __func__).raise();
    std::unique_ptr<SCOTCH_Strat, decltype(&SCOTCH_stratExit)> strategyCleanup(
      &strategy, &SCOTCH_stratExit);
    if (m_strategy && SCOTCH_stratGraphMap(&strategy, m_strategy->c_str()) != 0)
      Exception("SCOTCH_stratGraphMap failed", *this, __func__).raise();

    std::vector<SCOTCH_Num> partition(n, -1);
    if (SCOTCH_graphPart(&graph, static_cast<SCOTCH_Num>(numPartitions), &strategy,
          partition.data()) != 0)
      Exception("SCOTCH_graphPart failed", *this, __func__).raise();
    m_partition = std::move(partition);
    m_numPartitions = numPartitions;
  }

  size_t Partitioner::getPartition(Index index) const
  {
    if (index >= m_partition.size())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Partition assignment index " << index << " is outside assignment count "
        << m_partition.size() << '.' << Alert::Raise;
    }
    return m_partition[index];
  }
}
