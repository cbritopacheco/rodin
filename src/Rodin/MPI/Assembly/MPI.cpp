#include "MPI.h"

namespace Rodin::Assembly
{
  MPIIteration::MPIIteration(const Geometry::Mesh<Context::MPI>& mesh, Geometry::Region region)
    : m_mesh(mesh), m_region(region)
  {}

  Geometry::PolytopeIterator MPIIteration::getIterator() const
  {
    Geometry::PolytopeIterator it;
    switch (m_region)
    {
      case Geometry::Region::Cells:
        it = m_mesh.get().getCell();
        break;
      case Geometry::Region::Faces:
        it = m_mesh.get().getFace();
        break;
      case Geometry::Region::Boundary:
        it = m_mesh.get().getBoundary();
        break;
      case Geometry::Region::Interface:
        it = m_mesh.get().getInterface();
        break;
    }
    return it;
  }
}


