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
  /**
   * @brief Improves the mesh quality, mantaining the mean edge lenghts of the mesh.
   */
  class MeshOptimizer : public MMG5
  {
    public:
      MeshOptimizer() = default;

      /**
       * @brief Performs the optimization of the mesh.
       * @param[in, out] mesh Mesh to optimize
       * @note The mean of the edge lengths is preserved at the vertices.
       * Hence, if the edges passing through a vertex have very different
       * sizes, the resulting mesh may be very different from the initial one.
       */
      void optimize(Mesh<Traits::Serial>& mesh)
      {
        MMG5_pMesh mmgMesh = rodinToMesh(mesh);

        MMG5::setParameters(mmgMesh);

        bool isSurface = mesh.isSurface();
        switch (mmgMesh->dim)
        {
          case 2:
          {
            assert(!isSurface);
            optimizeMMG2D(mmgMesh);
            break;
          }
          case 3:
          {
            if (isSurface)
              optimizeMMGS(mmgMesh);
            else
              optimizeMMG3D(mmgMesh);
            break;
          }
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
      void optimizeMMG2D(MMG5_pMesh mesh);
      void optimizeMMG3D(MMG5_pMesh mesh);
      void optimizeMMGS(MMG5_pMesh mesh);
  };
}

#endif
