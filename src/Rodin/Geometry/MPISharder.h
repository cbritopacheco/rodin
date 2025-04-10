#ifndef RODIN_GEOMETRY_MPISHARDER_H
#define RODIN_GEOMETRY_MPISHARDER_H

#include <boost/mpi/config.hpp>

#include "Mesh.h"
#include "MPIMesh.h"

namespace Rodin::Geometry
{
  class MPISharder
  {
    public:
      MPISharder(const Context::MPI& context, LocalMesh& mesh);

      template <class CellPartition>
      Mesh<Context::MPI> shard(const CellPartition& partition);

    private:
      Context::MPI m_context;
      std::reference_wrapper<LocalMesh> m_mesh;
  };
}

#endif
