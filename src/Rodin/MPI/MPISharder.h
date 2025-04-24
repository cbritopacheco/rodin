#ifndef RODIN_GEOMETRY_MPISHARDER_H
#define RODIN_GEOMETRY_MPISHARDER_H

#include <boost/mpi/config.hpp>

#include "Rodin/Geometry/MeshPartitioner.h"

#include "MPIMesh.h"

namespace Rodin::Geometry
{
  /**
   * @class MPISharder
   * @brief Utility for distributing a global mesh across MPI ranks by
   *        splitting into per-rank shards, scattering them from a root,
   *        and gathering the local mesh on each rank.
   *
   * The typical usage is:
   * @code
   *   MPISharder sharder(ctx);
   *   auto mpiMesh = sharder.distribute(partitioner, rootRank);
   * @endcode
   */
  class MPISharder
  {
    public:
      /**
       * @brief Construct an MPISharder with the given MPI context.
       * @param context The MPI context (communicator and environment).
       */
      MPISharder(const Context::MPI& context)
        : m_context(context)
      {}

      /**
       * @brief One-step distribution: shard, scatter, and gather.
       *
       * Calls shard(), then scatter() on the root rank, then gather().
       *
       * @param p    The mesh partitioner defining per-cell owner ranks.
       * @param root The rank responsible for scattering shards.
       * @return The local MPI mesh built from the received shard.
       */
      Mesh<Context::MPI> distribute(Partitioner& p, int root)
      {
        return shard(p).scatter(root).gather(root);
      }

      /**
       * @brief Split the global mesh into per-rank shards (ghost layers included).
       *
       * Each cell of the global mesh is assigned to a shard builder
       * based on the partitioner, then each builder is finalized into a Shard.
       *
       * @param partitioner The mesh partitioner with `getCount()==comm.size()`.
       * @return Reference to this object for chaining.
       */
      MPISharder& shard(Partitioner& partitioner)
      {
        m_shards.clear();
        const auto& mesh = partitioner.getMesh();
        const size_t cellDim = mesh.getDimension();
        const size_t numShards = partitioner.getCount();
        assert(partitioner.getCount() == m_context.getCommunicator().size());
        assert(numShards > 0);
        std::vector<Shard::Builder> sbs(numShards);
        for (auto& sb : sbs)
          sb.initialize(mesh);
        for (auto it = mesh.getCell(); it; ++it)
        {
          const size_t partIdx = partitioner.getPartition(it->getIndex());
          sbs[partIdx].include(cellDim, it->getIndex());
        }
        m_shards.resize(numShards);
        for (size_t i = 0; i < numShards; i++)
          m_shards[i] = sbs[i].finalize();
        return *this;
      }

      /**
       * @brief Scatter shards from the root rank to all ranks.
       *
       * Only the root rank performs non-blocking sends of each shard
       * to its respective owner rank. Other ranks do nothing.
       *
       * @param root The rank that owns the global mesh and performs sends.
       * @return Reference to this object for chaining.
       */
      MPISharder& scatter(int root)
      {
        const auto& comm = m_context.getCommunicator();
        if (comm.rank() == root)
        {
          const int tag = m_context.getEnvironment().collectives_tag();
          std::vector<boost::mpi::request> requests;
          requests.reserve(m_shards.size());
          for (size_t i = 0; i < m_shards.size(); i++)
            requests.push_back(comm.isend(i, tag, m_shards[i]));
          for (auto& req : requests)
            req.wait();
        }
        return *this;
      }

      /**
       * @brief Gather the local shard on each rank.
       *
       * On the root rank, returns its own shard without MPI.
       * On other ranks, blocks in a matching MPI_Recv from the root.
       *
       * @param root The rank that originally scattered the shards.
       * @return The local MPI mesh built from the received shard.
       */
      Mesh<Context::MPI> gather(int root)
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

      Shard& getShard(size_t i)
      {
        assert(i < m_shards.size());
        return m_shards[i];
      }

      const Shard& getShard(size_t i) const
      {
        assert(i < m_shards.size());
        return m_shards[i];
      }

      const Context::MPI& getContext() const
      {
        return m_context;
      }

    private:
      Context::MPI m_context;
      std::vector<Shard> m_shards;
  };
}

#endif
