/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHOPTIMIZER_H
#define RODIN_EXTERNAL_MMG_MESHOPTIMIZER_H

#include "ForwardDecls.h"
#include "Mesh.h"
#include "MMG5.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Improves the mesh quality, mantaining the mean edge lenghts of the mesh.
   */
  class Optimize : public MMG5
  {
    public:
      Optimize() = default;

      /**
       * @brief Performs the optimization of the mesh.
       * @param[in, out] mesh Mesh to optimize
       * @note The mean of the edge lengths is preserved at the vertices.
       * Hence, if the edges passing through a vertex have very different
       * sizes, the resulting mesh may be very different from the initial one.
       */
      void optimize(MMG::Mesh& mesh)
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

      Optimize& setAngleDetection(bool b = true)
      {
        MMG5::setAngleDetection(b);
        return *this;
      }

      Optimize& setHMin(double hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      Optimize& setHMax(double hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      Optimize& setHausdorff(double hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      Optimize& setGradation(double hgrad)
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
