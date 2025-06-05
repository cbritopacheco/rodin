#include <boost/multi_array.hpp>

#include "Sharder.h"

namespace Rodin::Geometry
{
  MPISharder::Sharder(const Context::MPI& context)
    : m_context(context)
  {}

  MPIMesh MPISharder::distribute(Partitioner& p, int root)
  {
    return shard(p).scatter(root).gather(root);
  }

  MPISharder& Sharder<Context::MPI>::shard(Partitioner& partitioner)
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

    // Determine ownership of polytopes
    std::vector<std::vector<Boolean>> owned;
    owned.resize(cellDim + 1);
    for (size_t d = 0; d < cellDim + 1; d++)
      owned[d].resize(mesh.getPolytopeCount(d), false);
    for (Index i = 0; i < mesh.getCellCount(); i++)
    {
      const size_t partition = partitioner.getPartition(i);

      for (size_t d = 0; d < cellDim - 1; d++)
      {
        for (const Index idx : conn.getIncidence({ cellDim, d }, i))
        {
          if (owned[d][idx])
          {
            sbs[partition].include(d, idx, Shard::Flags::None);
          }
          else
          {
            sbs[partition].include(d, idx, Shard::Flags::Owned);
            owned[d][idx] = true;
          }
        }
      }

      assert(!owned[cellDim][i]);
      sbs[partition].include(cellDim, i, Shard::Flags::Owned);
      owned[cellDim][i] = true;
    }

    // Determine ghost polytopes
    for (Index i = 0; i < mesh.getCellCount(); i++)
    {
      const size_t partition = partitioner.getPartition(i);

      for (const Index nbr : conn.getIncidence({ cellDim, cellDim }, i))
      {
        if (partition == partitioner.getPartition(nbr))
          continue;

        for (size_t d = 0; d < cellDim - 1; d++)
        {
          for (const Index idx : conn.getIncidence({ cellDim, d }, nbr))
          {
            auto find = sbs[partition].getPolytopeMap(d).right.find(idx);
            if (find == sbs[partition].getPolytopeMap(d).right.end())
              sbs[partition].include(d, idx, Shard::Flags::Ghost);
          }
        }
        sbs[partition].include(cellDim, nbr, Shard::Flags::Ghost);
      }
    }

    m_shards.resize(numShards);
    for (size_t i = 0; i < numShards; i++)
      m_shards[i] = sbs[i].finalize();

    return *this;
  }

  MPISharder& MPISharder::scatter(int root)
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

  MPIMesh MPISharder::gather(int root)
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

  Shard& MPISharder::getShard(size_t i)
  {
    assert(i < m_shards.size());
    return m_shards[i];
  }

  const Shard& MPISharder::getShard(size_t i) const
  {
    assert(i < m_shards.size());
    return m_shards[i];
  }

  const Context::MPI& MPISharder::getContext() const
  {
    return m_context;
  }
}
