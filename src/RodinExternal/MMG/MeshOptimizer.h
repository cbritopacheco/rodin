/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHOPTIMIZER_H
#define RODIN_EXTERNAL_MMG_MESHOPTIMIZER_H

#include "Rodin/Mesh/Mesh.h"

#include "MMG5.h"

namespace Rodin::External::MMG
{
  class MeshOptimizer : public MMG5
  {
    public:
      MeshOptimizer() = default;

      void optimize(Mesh<Traits::Serial>& mesh)
      {
        MMG5_pMesh mmgMesh = rodinToMesh(mesh);

        MMG5::setParameters(mmgMesh);

        bool isSurface = mesh.isSurface();
        int retcode = MMG5_STRONGFAILURE;
        switch (mmgMesh->dim)
        {
          case 2:
          {
            assert(!isSurface);
            retcode = optimizeMMG2D(mmgMesh);
            break;
          }
          case 3:
          {
            if (isSurface)
              retcode = optimizeMMGS(mmgMesh);
            else
              retcode = optimizeMMG3D(mmgMesh);
            break;
          }
        }

        if (retcode != MMG5_SUCCESS)
        {
          Alert::Exception()
            << "Failed to optimize the mesh."
            << Alert::Raise;
        }

        mesh = meshToRodin(mmgMesh);
        destroyMesh(mmgMesh);
      }

      MeshOptimizer& setHMin(double hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      MeshOptimizer& setHMax(double hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      MeshOptimizer& setHausdorff(double hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      MeshOptimizer& setGradation(double hgrad)
      {
        MMG5::setGradation(hgrad);
        return *this;
      }

    private:
      int optimizeMMG2D(MMG5_pMesh mesh);
      int optimizeMMG3D(MMG5_pMesh mesh);
      int optimizeMMGS(MMG5_pMesh mesh);
  };
}

#endif
