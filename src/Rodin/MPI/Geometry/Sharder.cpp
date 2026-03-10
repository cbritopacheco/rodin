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

    const size_t D         = mesh.getDimension();
    const size_t numShards = partitioner.getCount();

    assert(numShards == static_cast<size_t>(comm.size()));
    assert(numShards > 0);

    const size_t numCells = mesh.getCellCount();

    // --------------------------------------------------------------------------
    // Cache global counts and cell partition once
    // --------------------------------------------------------------------------
    std::vector<size_t> polytopeCount(D + 1);
    for (size_t d = 0; d <= D; ++d)
      polytopeCount[d] = mesh.getPolytopeCount(d);

    std::vector<size_t> cellPartition(numCells);
#pragma omp parallel for schedule(static)
    for (Index c = 0; c < static_cast<Index>(numCells); ++c)
      cellPartition[c] = partitioner.getPartition(c);

    // Cache D->d incidences used in hot loops
    std::vector<const Incidence*> cellTo(D);
    for (size_t d = 0; d < D; ++d)
    {
      const auto& inc = conn.getIncidence(D, d);
      cellTo[d] = inc.empty() ? nullptr : &inc;
    }

    const auto& cellNbr = conn.getIncidence(D, D);

    // --------------------------------------------------------------------------
    // Builders
    // --------------------------------------------------------------------------
    std::vector<Shard::Builder> sbs(numShards);
    for (auto& sb : sbs)
      sb.initialize(mesh);

    // --------------------------------------------------------------------------
    // Global bookkeeping
    // owner[d][g] = owner shard of global polytope g in dimension d
    // local[d][g] = local index of g inside owner shard
    // halo[d][g]  = set of shards that need g as ghost
    // visited[d][g] = ownership assigned
    // --------------------------------------------------------------------------
    std::vector<std::vector<Index>> owner(D + 1);
    std::vector<std::vector<Index>> local(D + 1);
    std::vector<std::vector<IndexSet>> halo(D + 1);
    std::vector<std::vector<uint8_t>> visited(D + 1);

    for (size_t d = 0; d <= D; ++d)
    {
      owner[d].resize(polytopeCount[d]);
      local[d].resize(polytopeCount[d]);
      halo[d].resize(polytopeCount[d]);
      visited[d].assign(polytopeCount[d], uint8_t{0});
    }

    // Track only actually owned global entities so we do not rescan all globals.
    std::vector<std::vector<Index>> ownedGlobals(D + 1);
    for (size_t d = 0; d <= D; ++d)
      ownedGlobals[d].reserve(polytopeCount[d] / numShards + 1024);

    // --------------------------------------------------------------------------
    // Pass 1: assign owned entities from owned cells
    // --------------------------------------------------------------------------
    for (Index cellIdx = 0; cellIdx < static_cast<Index>(numCells); ++cellIdx)
    {
      const size_t part = cellPartition[cellIdx];
      auto& sb = sbs[part];

      for (size_t d = 0; d < D; ++d)
      {
        const Incidence* inc = cellTo[d];
        if (!inc)
          continue;

        for (const Index g : inc->at(cellIdx))
        {
          if (visited[d][g])
          {
            if (owner[d][g] != static_cast<Index>(part))
              halo[d][g].insert(part);

            const auto [localIdx, inserted] =
              sb.include({ d, g }, Shard::Flags::None);

            if (inserted)
              sb.getOwner(d)[localIdx] = owner[d][g];
          }
          else
          {
            const auto [localIdx, inserted] =
              sb.include({ d, g }, Shard::Flags::Owned);
            assert(inserted);

            owner[d][g]   = static_cast<Index>(part);
            local[d][g]   = localIdx;
            visited[d][g] = 1;
            ownedGlobals[d].push_back(g);
          }
        }
      }

      assert(!visited[D][cellIdx]);
      {
        const auto [localIdx, inserted] =
          sb.include({ D, cellIdx }, Shard::Flags::Owned);
        assert(inserted);

        owner[D][cellIdx]   = static_cast<Index>(part);
        local[D][cellIdx]   = localIdx;
        visited[D][cellIdx] = 1;
        ownedGlobals[D].push_back(cellIdx);
      }
    }

    // --------------------------------------------------------------------------
    // Pass 2: add ghost closure from neighboring cells in other partitions
    // --------------------------------------------------------------------------
    for (Index cellIdx = 0; cellIdx < static_cast<Index>(numCells); ++cellIdx)
    {
      const size_t part = cellPartition[cellIdx];
      auto& sb = sbs[part];

      for (const Index nbr : cellNbr.at(cellIdx))
      {
        const size_t nbrPart = cellPartition[nbr];
        if (nbrPart == part)
          continue;

        for (size_t d = 0; d < D; ++d)
        {
          const Incidence* inc = cellTo[d];
          if (!inc)
            continue;

          for (const Index g : inc->at(nbr))
          {
            auto find = sb.getPolytopeMap(d).right.find(g);
            if (find == sb.getPolytopeMap(d).right.end())
            {
              halo[d][g].insert(part);
              const auto [localIdx, inserted] =
                sb.include({ d, g }, Shard::Flags::Ghost);
              if (inserted)
                sb.getOwner(d)[localIdx] = owner[d][g];
            }
          }
        }

        halo[D][nbr].insert(part);
        const auto [localIdx, inserted] =
          sb.include({ D, nbr }, Shard::Flags::Ghost);
        if (inserted)
          sb.getOwner(D)[localIdx] = owner[D][nbr];
      }
    }

    // --------------------------------------------------------------------------
    // Move halo information only for actually owned entities
    // --------------------------------------------------------------------------
    for (size_t d = 0; d <= D; ++d)
    {
      for (const Index g : ownedGlobals[d])
      {
        auto& hm = sbs[owner[d][g]].getHalo(d);
        const Index k = local[d][g];

        auto [it, inserted] = hm.try_emplace(k, std::move(halo[d][g]));
        if (!inserted)
          it->second.insert(halo[d][g].begin(), halo[d][g].end());
      }
    }

    // --------------------------------------------------------------------------
    // Finalize shards in parallel
    // This is safe because each builder is independent.
    // --------------------------------------------------------------------------
    m_shards.resize(numShards);

#pragma omp parallel for schedule(dynamic, 1)
    for (Index p = 0; p < static_cast<Index>(numShards); ++p)
      m_shards[p] = sbs[p].finalize();

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
    boost::mpi::wait_all(reqs.begin(), reqs.end());
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
