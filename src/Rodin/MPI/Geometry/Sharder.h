#ifndef RODIN_GEOMETRY_MPISHARDER_H
#define RODIN_GEOMETRY_MPISHARDER_H

#include <boost/mpi/config.hpp>

#include "Rodin/Geometry/Sharder.h"
#include "Rodin/Geometry/MeshPartitioner.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
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
  template <>
  class Sharder<Context::MPI>
  {
    public:
      /**
       * @brief Construct an MPISharder with the given MPI context.
       * @param context The MPI context (communicator and environment).
       */
      Sharder(const Context::MPI& context);

      /**
       * @brief One-step distribution: shard, scatter, and gather.
       *
       * Calls shard(), then scatter() on the root rank, then gather().
       *
       * @param p    The mesh partitioner defining per-cell owner ranks.
       * @param root The rank responsible for scattering shards.
       * @return The local MPI mesh built from the received shard.
       */
      Mesh<Context::MPI> distribute(Partitioner& p, int root);

      /**
       * @brief Split the global mesh into per-rank shards (ghost layers included).
       *
       * Each cell of the global mesh is assigned to a shard builder
       * based on the partitioner, then each builder is finalized into a Shard.
       *
       * @param partitioner The mesh partitioner with `getCount()==comm.size()`.
       * @return Reference to this object for chaining.
       */
      Sharder& shard(Partitioner& partitioner);

      /**
       * @brief Scatter shards from the root rank to all ranks.
       *
       * Only the root rank performs non-blocking sends of each shard
       * to its respective owner rank. Other ranks do nothing.
       *
       * @param root The rank that owns the global mesh and performs sends.
       * @return Reference to this object for chaining.
       */
      Sharder& scatter(int root);

      /**
       * @brief Gather the local shard on each rank.
       *
       * On the root rank, returns its own shard without MPI.
       * On other ranks, blocks in a matching MPI_Recv from the root.
       *
       * @param root The rank that originally scattered the shards.
       * @return The local MPI mesh built from the received shard.
       */
      Mesh<Context::MPI> gather(int root);

      Shard& getShard(size_t i);

      const Shard& getShard(size_t i) const;

      const Context::MPI& getContext() const;

    private:
      Context::MPI m_context;
      std::vector<Shard> m_shards;
  };
}

namespace Rodin::MPI
{
  using Sharder = Geometry::Sharder<Context::MPI>;
}

#endif
