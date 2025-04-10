#include "MPISharder.h"

namespace Rodin::Geometry
{
  MPISharder::MPISharder(const Context::MPI& context, LocalMesh& mesh)
    : m_context(context),
      m_mesh(mesh)
  {}
}
