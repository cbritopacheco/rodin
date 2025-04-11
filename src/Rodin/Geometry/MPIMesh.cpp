#include "MPIMesh.h"

#ifdef RODIN_USE_MPI

namespace Rodin::Geometry
{
  MPIMesh::Builder::Builder(const Context::MPI& context)
    : m_context(context)
  {}

  MPIMesh::Mesh(const Context::MPI& context)
    : m_context(context)
  {}

  MPIMesh::shard(int rank, Shard&& shard)
  {
    const auto& comm = m_context.getCommunicator();
    if (comm.rank() == rank)
      comm.recv(rank, 0, m_shard);
    else
      comm.send(rank, 0, shard);
  }

  Index MPIMesh::getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId)
  {
    return 0;
  }
}

#endif
