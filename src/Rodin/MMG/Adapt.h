/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Adapt.h
 * @brief Metric-based mesh adaptation operator.
 */
#ifndef RODIN_EXTERNAL_MMG_ADAPT_H
#define RODIN_EXTERNAL_MMG_ADAPT_H

#include "ForwardDecls.h"

#include "Mesh.h"

#include "MMG5.h"

namespace Rodin::MMG
{
  /**
   * @brief Performs anisotropic/isotropic remeshing from a nodal size map.
   *
   * This class wraps MMG adaptation routines and dispatches automatically to:
   * - MMG2D for planar 2D meshes,
   * - MMG3D for volumetric 3D meshes,
   * - MMGS for surface meshes embedded in 3D.
   *
   * The scalar size-map field encodes target local edge length and is sampled
   * at mesh vertices.
   *
   * Parameters inherited from @ref MMG5 (hmin/hmax/hausdorff/gradation/angle)
   * can be configured before calling @ref adapt.
   */
  class Adapt : public MMG5
  {
    public:
      /**
       * @brief Constructs an adaptation operator with default MMG settings.
       */
      Adapt() = default;

      /**
       * @brief Adapts a mesh according to a scalar nodal size map.
       * @param[in, out] mesh Mesh to remesh in place.
       * @param[in] sizeMap Scalar P1 field prescribing local target sizes.
       *
       * The input mesh is converted to MMG native structures, remeshed, then
       * converted back to @ref MMG::Mesh.
       */
      void adapt(MMG::Mesh& mesh, const RealGridFunction& sizeMap);

      /**
       * @brief Enables/disables ridge angle detection.
       * @param[in] b If `true`, enable detection.
       * @returns Reference to this object.
       */
      Adapt& setAngleDetection(bool b = true)
      {
        MMG5::setAngleDetection(b);
        return *this;
      }

      /**
       * @brief Sets the minimum edge size constraint.
       * @param[in] hmin Minimum size.
       * @returns Reference to this object.
       */
      Adapt& setHMin(double hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      /**
       * @brief Sets the maximum edge size constraint.
       * @param[in] hmax Maximum size.
       * @returns Reference to this object.
       */
      Adapt& setHMax(double hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      /**
       * @brief Sets the Hausdorff tolerance controlling boundary approximation.
       * @param[in] hausd Hausdorff value.
       * @returns Reference to this object.
       */
      Adapt& setHausdorff(double hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      /**
       * @brief Sets the allowed gradation between adjacent edge sizes.
       * @param[in] hgrad Gradation ratio.
       * @returns Reference to this object.
       */
      Adapt& setGradation(double hgrad)
      {
        MMG5::setGradation(hgrad);
        return *this;
      }

    private:
      /**
       * @brief Runs MMG2D adaptation kernel.
       * @param[in, out] mesh MMG2D mesh.
       * @param[in, out] size MMG size-map solution.
       * @returns MMG return code.
       */
      ReturnCode adaptMMG2D(MMG5_pMesh mesh, MMG5_pSol size);
      /**
       * @brief Runs MMG3D adaptation kernel.
       * @param[in, out] mesh MMG3D mesh.
       * @param[in, out] size MMG size-map solution.
       * @returns MMG return code.
       */
      ReturnCode adaptMMG3D(MMG5_pMesh mesh, MMG5_pSol size);
      /**
       * @brief Runs MMGS adaptation kernel for surface meshes.
       * @param[in, out] mesh MMGS mesh.
       * @param[in, out] size MMG size-map solution.
       * @returns MMG return code.
       */
      ReturnCode adaptMMGS(MMG5_pMesh mesh, MMG5_pSol size);
  };
}

#endif
