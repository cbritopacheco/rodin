/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MeshOptimizer.h
 * @brief Mesh-quality optimization operator.
 */
#ifndef RODIN_EXTERNAL_MMG_MESHOPTIMIZER_H
#define RODIN_EXTERNAL_MMG_MESHOPTIMIZER_H

#include "Rodin/Alert/MemberFunctionWarning.h"

#include "ForwardDecls.h"
#include "Mesh.h"
#include "MMG5.h"

namespace Rodin::MMG
{
  /**
   * @brief Improves mesh quality while preserving characteristic local sizes.
   *
   * This operator runs MMG optimization mode without prescribing an external
   * metric field. It is useful for improving element shapes while keeping the
   * mesh close to the current local edge-size distribution.
   *
   * Typical usage is post-processing after adaptation/imported meshes to remove
   * poorly shaped elements while preserving geometric features.
   */
  class Optimizer : public MMG5
  {
    public:
      /**
       * @brief Default constructor.
       */
      Optimizer() = default;

      /**
       * @brief Optimizes mesh quality in place.
       * @param[in, out] mesh Mesh to optimize.
       *
       * Dispatches to MMG2D/MMG3D/MMGS according to mesh type.
       *
       * @note The mean edge length around vertices is preserved. If local sizes
       * vary strongly around a vertex, the optimized mesh can still differ
       * significantly from the initial one.
       */
      void optimize(MMG::Mesh& mesh)
      {
        if (mesh.isEmpty())
        {
          Alert::MemberFunctionWarning(*this, __func__)
            << "The MMG::Mesh object is empty. No action performed."
            << Alert::Raise;
          return;
        }

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
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to optimize the mesh."
            << Alert::Raise;
        }

        mesh = meshToRodin(mmgMesh);
        destroyMesh(mmgMesh);
      }

      /**
       * @brief Enables/disables ridge angle detection.
       * @param[in] b If `true`, enable detection.
       * @returns Reference to this object.
       */
      Optimizer& setAngleDetection(bool b = true)
      {
        MMG5::setAngleDetection(b);
        return *this;
      }

      /**
       * @brief Sets the minimum edge size constraint.
       * @param[in] hmin Minimum size.
       * @returns Reference to this object.
       */
      Optimizer& setHMin(double hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      /**
       * @brief Sets the maximum edge size constraint.
       * @param[in] hmax Maximum size.
       * @returns Reference to this object.
       */
      Optimizer& setHMax(double hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      /**
       * @brief Sets the Hausdorff tolerance controlling boundary approximation.
       * @param[in] hausd Hausdorff value.
       * @returns Reference to this object.
       */
      Optimizer& setHausdorff(double hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      /**
       * @brief Sets the allowed gradation between adjacent edge sizes.
       * @param[in] hgrad Gradation ratio.
       * @returns Reference to this object.
       */
      Optimizer& setGradation(double hgrad)
      {
        MMG5::setGradation(hgrad);
        return *this;
      }

    private:
      /**
       * @brief Runs MMG2D optimization kernel.
       * @param[in, out] mesh MMG2D mesh.
       * @returns MMG return code.
       */
      int optimizeMMG2D(MMG5_pMesh mesh);
      /**
       * @brief Runs MMG3D optimization kernel.
       * @param[in, out] mesh MMG3D mesh.
       * @returns MMG return code.
       */
      int optimizeMMG3D(MMG5_pMesh mesh);
      /**
       * @brief Runs MMGS optimization kernel for surface meshes.
       * @param[in, out] mesh MMGS mesh.
       * @returns MMG return code.
       */
      int optimizeMMGS(MMG5_pMesh mesh);
  };
}

#endif
