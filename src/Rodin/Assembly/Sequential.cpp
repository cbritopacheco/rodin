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

  Geometry::PolytopeIterator
  SequentialIteration<Geometry::Mesh<Context::Local>>::getIterator() const
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

  Geometry::Polytope SequentialIteration<Geometry::Mesh<Context::Local>>::getPolytope(
    Index i) const
  {
    return Geometry::Polytope(getDimension(), i, m_mesh.get());
  }

  size_t SequentialIteration<Geometry::Mesh<Context::Local>>::getDimension() const
  {
    switch (m_region)
    {
      case Geometry::Region::Cells:
        return m_mesh.get().getDimension();
      case Geometry::Region::Boundary:
      case Geometry::Region::Interface:
      case Geometry::Region::Faces:
        return m_mesh.get().getDimension() - 1;
    }
    assert(false);
    return 0;
  }

  size_t SequentialIteration<Geometry::Mesh<Context::Local>>::getCount() const
  {
    return m_region == Geometry::Region::Cells ? m_mesh.get().getCellCount()
                                               : m_mesh.get().getFaceCount();
  }

  bool SequentialIteration<Geometry::Mesh<Context::Local>>::filter(Index i) const
  {
    switch (m_region)
    {
      case Geometry::Region::Cells:
      case Geometry::Region::Faces:
        return true;
      case Geometry::Region::Boundary:
        return m_mesh.get().isBoundary(i);
      case Geometry::Region::Interface:
        return m_mesh.get().isInterface(i);
    }
    assert(false);
    return false;
  }
}
