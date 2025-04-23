/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>

#include "Multithreaded.h"

namespace Rodin::Assembly
{
  namespace Internal
  {
    MultithreadedIteration::MultithreadedIteration(const Geometry::MeshBase& mesh, Variational::Integrator::Region region)
      : m_mesh(mesh), m_region(region)
    {}

    Geometry::PolytopeIterator MultithreadedIteration::getIterator(Index i) const
    {
      Geometry::PolytopeIterator it;
      switch (m_region)
      {
        case Variational::Integrator::Region::Cells:
        {
          it = m_mesh.get().getCell(i);
          return it;
        }
        case Variational::Integrator::Region::Faces:
        case Variational::Integrator::Region::Boundary:
        case Variational::Integrator::Region::Interface:
        {
          it = m_mesh.get().getFace(i);
          return it;
        }
      }
      assert(false);
      return it;
    }

    size_t MultithreadedIteration::getDimension() const
    {
      switch (m_region)
      {
        case Variational::Integrator::Region::Cells:
        {
          return m_mesh.get().getDimension();
        }
        case Variational::Integrator::Region::Faces:
        case Variational::Integrator::Region::Boundary:
        case Variational::Integrator::Region::Interface:
        {
          return m_mesh.get().getDimension() - 1;
        }
      }
      assert(false);
      return 0;
    }

    size_t MultithreadedIteration::getCount() const
    {
      switch (m_region)
      {
        case Variational::Integrator::Region::Cells:
        {
          return m_mesh.get().getCellCount();
        }
        case Variational::Integrator::Region::Faces:
        case Variational::Integrator::Region::Boundary:
        case Variational::Integrator::Region::Interface:
        {
          return m_mesh.get().getFaceCount();
        }
      }
      assert(false);
      return 0;
    }

    bool MultithreadedIteration::filter(Index i) const
    {
      switch (m_region)
      {
        case Variational::Integrator::Region::Faces:
        case Variational::Integrator::Region::Cells:
        {
          return true;
        }
        case Variational::Integrator::Region::Boundary:
        {
          return m_mesh.get().isBoundary(i);
        }
        case Variational::Integrator::Region::Interface:
        {
          return m_mesh.get().isInterface(i);
        }
      }
      assert(false);
      return false;
    }
  }
}

