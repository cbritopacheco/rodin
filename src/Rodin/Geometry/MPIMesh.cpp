#include "MPIMesh.h"

#ifdef RODIN_USE_MPI

namespace Rodin::Geometry
{
  MPIMesh::Builder& MPIMesh::Builder::initialize(const Context::MPI& context, Shard&& shard)
  {
    m_context = std::move(context);
    m_shard = std::move(shard);
    return *this;
  }

  // MPIMesh MPIMesh::Builder::finalize()
  // {
  //   MPIMesh mesh(m_context);
  //   mesh.m_shard = std::move(m_shard);
  //   return mesh;
  // }

  Index MPIMesh::getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId)
  {
    return 0;
  }
}

#endif
