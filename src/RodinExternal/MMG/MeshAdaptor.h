/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHADAPTOR_H
#define RODIN_EXTERNAL_MMG_MESHADAPTOR_H

#include "Rodin/Mesh/Mesh.h"
#include "Rodin/Variational/GridFunction.h"

#include "MMG5.h"

namespace Rodin::External::MMG
{
  class MeshAdaptor : public MMG5
  {
    public:
      MeshAdaptor() = default;

      /**
       * @brief Performs the optimization of the mesh.
       * @param[in, out] mesh Mesh to optimize
       * @note The mean of the edge lengths is preserved at the vertices.
       * Hence, if the edges passing through a vertex have very different
       * sizes, the resulting mesh may be very different from the initial one.
       */
      template <class FES>
      void adapt(Mesh<Context::Serial>& mesh, const Variational::GridFunction<FES>& gf)
      {
        MMG5_pMesh mmgMesh = rodinToMesh(mesh);

        MMG5::setParameters(mmgMesh);

        assert(gf.getFiniteElementSpace().getVectorDimension() == 1);
        MMG5_pSol sol = createSolution(mmgMesh, 1);
        copySolution(gf, sol);

        bool isSurface = mesh.isSurface();
        int retcode = MMG5_STRONGFAILURE;
        switch (mmgMesh->dim)
        {
          case 2:
          {
            assert(!isSurface);
            retcode = adaptMMG2D(mmgMesh, sol);
            break;
          }
          case 3:
          {
            if (isSurface)
              retcode = adaptMMGS(mmgMesh, sol);
            else
              retcode = adaptMMG3D(mmgMesh, sol);
            break;
          }
        }

        if (retcode != MMG5_SUCCESS)
        {
          Alert::Exception()
            << "Failed to adapt the mesh."
            << Alert::Raise;
        }

        mesh = meshToRodin(mmgMesh);

        destroySolution(sol);
        destroyMesh(mmgMesh);
      }

      MeshAdaptor& setHMin(double hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      MeshAdaptor& setHMax(double hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      MeshAdaptor& setHausdorff(double hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      MeshAdaptor& setGradation(double hgrad)
      {
        MMG5::setGradation(hgrad);
        return *this;
      }

    private:
      int adaptMMG2D(MMG5_pMesh mesh, MMG5_pSol sol)
      {
        return 0;
      }

      int adaptMMG3D(MMG5_pMesh mesh, MMG5_pSol sol);

      int adaptMMGS(MMG5_pMesh mesh, MMG5_pSol sol)
      {
        return 0;
      }
  };
}

#endif
