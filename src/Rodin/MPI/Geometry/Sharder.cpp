#include <cassert>
#include <boost/multi_array.hpp>

#include "Rodin/Serialization/BitSet.h"
#include "Rodin/Serialization/FlatMap.h"
#include "Rodin/Serialization/FlatSet.h"
#include "Rodin/Serialization/Optional.h"

#include "Sharder.h"

namespace Rodin::Geometry
{
  Sharder<Context::MPI>::Sharder(const Context::MPI& context)
    : m_context(context)
  {}

  MPIMesh Sharder<Context::MPI>::distribute(Partitioner& p, int root)
  {
    return shard(p).scatter(root).gather(root);
  }

  Sharder<Context::MPI>& Sharder<Context::MPI>::shard(Partitioner& partitioner)
  {
    m_shards.clear();
    const auto& mesh = partitioner.getMesh();
    const auto& conn = mesh.getConnectivity();
    const size_t cellDim = mesh.getDimension();
    const size_t numShards = partitioner.getCount();
    assert(m_context.getCommunicator().size() >= 0);
    assert(partitioner.getCount() == static_cast<size_t>(m_context.getCommunicator().size()));
    assert(numShards > 0);
    std::vector<Shard::Builder> sbs(numShards);
    for (auto& sb : sbs)
      sb.initialize(mesh);
    std::vector<std::vector<Index>> owner(cellDim + 1);     // Owner of each polytope
    std::vector<std::vector<Index>> local(cellDim + 1);     // Owning shard index for each polytope
    std::vector<std::vector<IndexSet>> halo(cellDim + 1);   // Halo information for each polytope
    std::vector<std::vector<uint8_t>> visited(cellDim + 1); // Visited polytopes

    for (size_t d = 0; d < cellDim + 1; d++)
    {
      owner[d].resize(mesh.getPolytopeCount(d));
      local[d].resize(mesh.getPolytopeCount(d));
      halo[d].resize(mesh.getPolytopeCount(d));
      visited[d].resize(mesh.getPolytopeCount(d), false);
    }

    for (Index cellIdx = 0; cellIdx < mesh.getCellCount(); cellIdx++)
    {
      const size_t partition = partitioner.getPartition(cellIdx);
      for (size_t d = 0; d <= cellDim - 1; d++)
      {
        const auto& inc = conn.getIncidence(cellDim, d);
        if (inc.empty())
          continue;
        for (const Index& idx : inc.at(cellIdx))
        {
          if (visited[d][idx])
          {
            if (owner[d][idx] != partition)
              halo[d][idx].insert(partition);
            const auto [localIdx, inserted] =
              sbs[partition].include({ d, idx }, Shard::Flags::None);
            if (inserted)
              sbs[partition].getOwner(d)[localIdx] = owner[d][idx];
          }
          else
          {
            const auto [localIdx, inserted] =
              sbs[partition].include({ d, idx }, Shard::Flags::Owned);
            assert(inserted);
            owner[d][idx] = partition;
            local[d][idx] = localIdx;
            visited[d][idx] = true;
          }
        }
      }
      assert(!visited[cellDim][cellIdx]);
      const auto [localIdx, inserted] =
        sbs[partition].include({ cellDim, cellIdx }, Shard::Flags::Owned);
      assert(inserted);
      owner[cellDim][cellIdx] = partition;
      local[cellDim][cellIdx] = localIdx;
      visited[cellDim][cellIdx] = true;
    }

    for (Index cellIdx = 0; cellIdx < mesh.getCellCount(); cellIdx++)
    {
      const size_t partition = partitioner.getPartition(cellIdx);
      for (const Index& nbr : conn.getIncidence({ cellDim, cellDim }, cellIdx))
      {
        if (partition == partitioner.getPartition(nbr))
          continue;
        for (size_t d = 0; d <= cellDim - 1; d++)
        {
          const auto& inc = conn.getIncidence(cellDim, d);
          if (inc.empty())
            continue;
          for (const Index idx : inc.at(nbr))
          {
            auto find = sbs[partition].getPolytopeMap(d).right.find(idx);
            if (find == sbs[partition].getPolytopeMap(d).right.end())
            {
              halo[d][idx].insert(partition);
              const auto [localIdx, inserted] =
                sbs[partition].include({ d, idx }, Shard::Flags::Ghost);
              if (inserted)
                sbs[partition].getOwner(d)[localIdx] = owner[d][idx];
            }
          }
        }
        halo[cellDim][nbr].insert(partition);
        const auto [localIdx, inserted] =
          sbs[partition].include({ cellDim, nbr }, Shard::Flags::Ghost);
        if (inserted)
          sbs[partition].getOwner(cellDim)[localIdx] = owner[cellDim][nbr];
      }
    }

    // Now move the halo information into the shards
    for (size_t d = 0; d < cellDim + 1; d++)
    {
      for (size_t i = 0; i < mesh.getPolytopeCount(d); i++)
      {
        if (!visited[d][i])
          continue;
        auto& hm = sbs[owner[d][i]].getHalo(d);
        const size_t k = local[d][i];
        auto [it, inserted] = hm.try_emplace(k, std::move(halo[d][i]));
        if (!inserted)
          it->second.insert(halo[d][i].begin(), halo[d][i].end());
      }
    }

    m_shards.resize(numShards);
    for (size_t i = 0; i < numShards; i++)
      m_shards[i] = sbs[i].finalize();

    return *this;
  }

  Sharder<Context::MPI>& Sharder<Context::MPI>::scatter(int root)
  {
    const auto& comm = m_context.getCommunicator();
    const int tag = m_context.getEnvironment().collectives_tag();
    std::vector<boost::mpi::request> reqs(m_shards.size());
    for (size_t i = 0; i < m_shards.size(); i++)
    {
      assert(root >= 0);
      if (i == static_cast<size_t>(root))
        continue;
      reqs[i] = comm.isend(i, tag, m_shards[i]);
    }
    for (auto& req : reqs)
      req.wait();
    return *this;
  }

  MPIMesh Sharder<Context::MPI>::gather(int root)
  {
    const auto& comm = m_context.getCommunicator();
    const int tag = m_context.getEnvironment().collectives_tag();
    if (comm.rank() == root)
    {
      return Mesh<Context::MPI>::Builder(m_context).initialize(std::move(m_shards[root]))
                                                   .finalize();
    }
    else
    {
      Shard s;
      comm.recv(root, tag, s);
      Mesh<Context::MPI>::Builder build(m_context);
      build.initialize(std::move(s));
      return build.finalize();
    }
  }

  Shard& Sharder<Context::MPI>::getShard(size_t i)
  {
    assert(i < m_shards.size());
    return m_shards[i];
  }

  const Shard& Sharder<Context::MPI>::getShard(size_t i) const
  {
    assert(i < m_shards.size());
    return m_shards[i];
  }

  const Context::MPI& Sharder<Context::MPI>::getContext() const
  {
    return m_context;
  }
}
