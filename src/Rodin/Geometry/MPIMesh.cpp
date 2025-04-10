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

  Index MPIMesh::getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId)
  {
    return 0;
  }
}

#endif
