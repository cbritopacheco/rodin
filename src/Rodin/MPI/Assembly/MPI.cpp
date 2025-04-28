#include "MPI.h"

namespace Rodin::Assembly::Internal
{
  MPIIteration::MPIIteration(const Geometry::Mesh<Context::MPI>& mesh, Variational::Integrator::Region region)
    : m_mesh(mesh), m_region(region)
  {}

  Geometry::PolytopeIterator MPIIteration::getIterator() const
  {
    Geometry::PolytopeIterator it;
    switch (m_region)
    {
      case Variational::Integrator::Region::Cells:
      {
        it = m_mesh.get().getShard().getCell();
        break;
      }
      case Variational::Integrator::Region::Faces:
      {
        it = m_mesh.get().getShard().getFace();
        break;
      }
      case Variational::Integrator::Region::Boundary:
      {
        it = m_mesh.get().getShard().getBoundary();
        break;
      }
      case Variational::Integrator::Region::Interface:
      {
        it = m_mesh.get().getShard().getInterface();
        break;
      }
    }
    return it;
  }
}


