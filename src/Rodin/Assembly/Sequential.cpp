/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/Mesh.h"

#include "Sequential.h"

namespace Rodin::Assembly
{
  SequentialIteration<Geometry::Mesh<Context::Local>>
  ::SequentialIteration(const Geometry::Mesh<Context::Local>& mesh, const Geometry::Region& region)
    : m_mesh(mesh), m_region(region)
  {}

  Geometry::PolytopeIterator SequentialIteration<Geometry::Mesh<Context::Local>>::getIterator() const
  {
    Geometry::PolytopeIterator it;
    switch (m_region)
    {
      case Geometry::Region::Cells:
      {
        it = m_mesh.get().getCell();
        break;
      }
      case Geometry::Region::Faces:
      {
        it = m_mesh.get().getFace();
        break;
      }
      case Geometry::Region::Boundary:
      {
        it = m_mesh.get().getBoundary();
        break;
      }
      case Geometry::Region::Interface:
      {
        it = m_mesh.get().getInterface();
        break;
      }
    }
    return it;
  }
}


