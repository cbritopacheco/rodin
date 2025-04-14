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
      MPISharder(const Context::MPI& context)
        : m_context(context),
          m_root(0),
          m_tag(0)
      {}

      void shard(Partitioner& partitioner)
      {
        const auto& comm = m_context.getCommunicator();
        const auto& mesh = partitioner.getMesh();
        const size_t cellDim = mesh.getDimension();
        const size_t numShards = partitioner.getCount();

        Mesh<Context::MPI>::Builder build;
        if (comm.rank() == m_root)
        {
          std::vector<Shard::Builder> sbs(numShards);
          for (auto& sb : sbs)
            sb.initialize(mesh);

          for (auto it = mesh.getCell(); it; ++it)
          {
            const size_t partIdx = partitioner.getPartition(it->getIndex());
            sbs[partIdx].include(cellDim, it->getIndex());
          }

          std::vector<Shard> shards(numShards);
          std::vector<boost::mpi::request> requests;
          requests.reserve(numShards - 1);
          for (size_t i = 0; i < shards.size(); i++)
          {
            assert(m_root >= 0);
            if (i == static_cast<size_t>(m_root))
              shards[i] = sbs[i].finalize();
            else
              requests.push_back(comm.isend(i, m_tag, sbs[i].finalize()));
          }
          for (auto& req : requests)
            req.wait();
          build.initialize(m_context, std::move(shards[m_root]));
        }
        else
        {
          Shard s;
          boost::mpi::request request = comm.irecv(m_root, m_tag, s);
          request.wait();
          Mesh<Context::MPI>::Builder build;
          build.initialize(m_context, std::move(s));
        }
        // return build.finalize();
      }

    private:
      Context::MPI m_context;
      int m_root;
      int m_tag;
  };
}

#endif
