/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_ADAPT_H
#define RODIN_EXTERNAL_MMG_ADAPT_H

#include "ForwardDecls.h"

#include "Mesh.h"

#include "MMG5.h"

namespace Rodin::External::MMG
{
  class Adapt : public MMG5
  {
    public:
      Adapt() = default;

      void adapt(MMG::Mesh& mesh, const ScalarGridFunction& sizeMap);

      Adapt& setAngleDetection(bool b = true)
      {
        MMG5::setAngleDetection(b);
        return *this;
      }

      Adapt& setHMin(double hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      Adapt& setHMax(double hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      Adapt& setHausdorff(double hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      Adapt& setGradation(double hgrad)
      {
        MMG5::setGradation(hgrad);
        return *this;
      }

    private:
      ReturnCode adaptMMG2D(MMG5_pMesh mesh, MMG5_pSol size);
      ReturnCode adaptMMG3D(MMG5_pMesh mesh, MMG5_pSol size);
      ReturnCode adaptMMGS(MMG5_pMesh mesh, MMG5_pSol size);
  };
}

#endif

