/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Sequential.h"

namespace Rodin::Assembly
{
  namespace Internal
  {
    Geometry::PolytopeIterator getRegionIterator(
        const Geometry::MeshBase& mesh, Variational::Integrator::Region region)
    {
      Geometry::PolytopeIterator it;
      switch (region)
      {
        case Variational::Integrator::Region::Cells:
        {
          it = mesh.getCell();
          break;
        }
        case Variational::Integrator::Region::Faces:
        {
          it = mesh.getFace();
          break;
        }
        case Variational::Integrator::Region::Boundary:
        {
          it = mesh.getBoundary();
          break;
        }
        case Variational::Integrator::Region::Interface:
        {
          it = mesh.getInterface();
          break;
        }
      }
      return it;
    }
  }
}


