#ifndef RODIN_GEOMETRY_MPISHARDER_H
#define RODIN_GEOMETRY_MPISHARDER_H

#include <boost/mpi/config.hpp>

#include "Mesh.h"
#include "MPIMesh.h"
#include "MeshPartitioner.h"

namespace Rodin::Geometry
{
  class MPISharder
  {
    public:
      MPISharder(const Context::MPI& context);

      Mesh<Context::MPI> shard(Partitioner& partitioner)
      {
        const auto& localMesh = partitioner.getMesh();
        const size_t cellDim = localMesh.getDimension();
        const size_t numShards = partitioner.getCount();

        Mesh<Context::MPI>::Builder build;
        build.initialize(localMesh.getSpaceDimension(), m_context);

        std::vector<Shard::Builder> shards;
        shards.resize(numShards);

        for (auto& shard : shards)
          shard.initialize(localMesh);

        for (auto it = localMesh.getCell(); it; ++it)
        {
          const size_t partIdx = partitioner.getPartition(it->getIndex());
          shards[partIdx].include(cellDim, it->getIndex());
        }

        // for (const auto& shard : shards)
        // {
        //   const auto& l = shard.finalize();
        //   const size_t numCells = localShard.getCellCount();
        //   for (size_t i = 0; i < numCells; i++)
        //   {
        //     const auto& polytope = localShard.getPolytope(cellDim, i);
        //     const auto& globalIdx = polytope.getGlobalIndex();
        //     build.polytope(polytope.getGeometryType(), globalIdx);
        //   }
        // }

      }

    private:
      Context::MPI m_context;
      std::reference_wrapper<LocalMesh> m_mesh;
  };
}

#endif
