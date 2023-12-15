/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Sequential.h"

namespace Rodin::Assembly::Internal
{
  SequentialIteration::SequentialIteration(const Geometry::MeshBase& mesh, Variational::Integrator::Region region)
    : m_mesh(mesh), m_region(region)
  {}

  Geometry::PolytopeIterator SequentialIteration::getIterator() const
  {
    Geometry::PolytopeIterator it;
    switch (m_region)
    {
      case Variational::Integrator::Region::Cells:
      {
        it = m_mesh.get().getCell();
        break;
      }
      case Variational::Integrator::Region::Faces:
      {
        it = m_mesh.get().getFace();
        break;
      }
      case Variational::Integrator::Region::Boundary:
      {
        it = m_mesh.get().getBoundary();
        break;
      }
      case Variational::Integrator::Region::Interface:
      {
        it = m_mesh.get().getInterface();
        break;
      }
    }
    return it;
  }
}


