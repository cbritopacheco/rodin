#ifndef RODIN_MPI_GEOMETRY_SHARDER_H
#define RODIN_MPI_GEOMETRY_SHARDER_H

/**
 * @file
 * @brief Mesh sharding and distribution utilities for MPI contexts.
 */

#include <boost/mpi/config.hpp>

#include "Rodin/Geometry/Sharder.h"
#include "Rodin/Geometry/MeshPartitioner.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Utility for distributing a global mesh across MPI ranks by
   * splitting into per-rank shards, scattering them from a root, and gathering
   * the local mesh on each rank.
   */
  template <>
  class Sharder<Context::MPI> : public SharderBase<Context::MPI>
  {
    public:
      /**
       * @brief Base sharder interface specialized on MPI context.
       */
      using Parent = SharderBase<Context::MPI>;

      /**
       * @brief Construct a Sharder with the given MPI context.
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
  };
}

namespace Rodin::MPI
{
  /**
   * @brief Convenience alias for the MPI mesh sharder specialization.
   */
  using Sharder = Geometry::Sharder<Context::MPI>;
}

#endif
