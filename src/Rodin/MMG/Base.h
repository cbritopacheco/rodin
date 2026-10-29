/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Base.h
 * @brief Shared configuration for MMG-backed remeshing operations.
 */
#ifndef RODIN_MMG_BASE_H
#define RODIN_MMG_BASE_H

#include <boost/filesystem.hpp>

#include "Rodin/Math.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Utility/IsSpecialization.h"

#include "Common.h"
#include "Configure.h"
#include "ForwardDecls.h"

namespace Rodin::MMG
{
  /**
   * @brief Exception type for MMG wrapper failures.
   * @tparam FuncName Type of the function-name object passed by caller.
   *
   * This exception prepends a standardized "class/function" prefix to MMG
   * diagnostic messages, then can be streamed with additional details before
   * raising through Rodin's alert system.
   */
  template <class FuncName>
  class MMGException : public Alert::Exception
  {
    public:
      /**
       * @brief Builds an MMG wrapper exception message prefix.
       * @param[in] funcName Function identifier used in the diagnostic message.
       */
      MMGException(const FuncName& funcName)
      {
        const auto& className = boost::typeindex::type_id_with_cvr<Base>().pretty_name();
        *this << "In " << Alert::Identifier::Function(funcName) << " of class "
              << Alert::Identifier::Class(className) << ": ";
      }
  };

  /**
   * @brief Return-code type used by MMG C API calls.
   *
   * MMG kernels conventionally return `MMG5_SUCCESS`, `MMG5_LOWFAILURE`, or
   * `MMG5_STRONGFAILURE`.
   */
  using ReturnCode = int;

  /**
   * @brief Shared configuration for MMG-backed remeshing operations.
   *
   * High-level operators such as @ref Adapt, @ref Optimizer, and
   * @ref LevelSetDiscretizer inherit the common MMG controls.
   */
  class Base
  {
    public:
      /// @brief Constructs the shared MMG remeshing configuration.
      Base();

      /**
       * @brief Enables/disables ridge angle detection.
       * @param[in] b If `true`, enable angle detection.
       * @returns Reference to this object.
       */
      Base& setAngleDetection(bool b = true)
      {
        m_ridgeDetection = b;
        return *this;
      }

      /**
       * @brief Sets the minimal edge size.
       *
       * @param[in] hmin Minimal edge size.
       *
       * @returns Reference to self (for method chaining)
       *
       * The `hmin` option allows to truncate the edge sizes to be lower than the
       * `hmax` one.
       *
       * The default values for this parameters are computed from the mesh
       * bounding box or, if provided, from the given metric.
       *
       * - Without metric, the minimal edge size is set to 0.01 of the bounding
       * box size.
       *
       * - With metric, the minimal edge size is set to 0.1 of the
       * smallest prescribed size.
       *
       * @see setHMax(Real)
       */
      Base& setHMin(Real hmin)
      {
        m_hmin = hmin;
        return *this;
      }

      /**
       * @brief Sets the maximal edge size parameter.
       *
       * @param[in] hmax Maximal edge size.
       *
       * @returns Reference to self (for method chaining)
       *
       * The `hmax` option allows to truncate the edge sizes to be greater than
       * the `hmin` parameter.
       *
       * The default values for this parameters are computed from the mesh
       * bounding box or, if provided, from the given metric.
       *
       * - Without metric, the maximal edge size is set to two times the
       *  bounding box size.
       *
       * - With metric, the maximal one is set to 10 times the maximal
       *  prescribed size.
       *
       * @see setHMin(Real)
       */
      Base& setHMax(Real hmax)
      {
        m_hmax = hmax;
        return *this;
      }

      /**
       * @brief Sets the Hausdorff parameter.
       *
       * @param[in] hausd Hausdorff parameter.
       *
       * @returns Reference to self (for method chaining)
       *
       * The Hausdorff parameter controls the boundary approximation.  It
       * imposes the maximal distance between the piecewise linear
       * representation of the boundary and the reconstructed ideal boundary.
       * Thus, a low Hausdorff parameter leads to the refinement of high
       * curvature areas.
       *
       * By default, the Hausdorff value is set to 0.01, which is a suitable
       * value for an object of size 1 in each direction. For smaller (resp.
       * larger) objects, you may need to decrease (resp. increase) the
       * Hausdorff parameter.
       *
       */
      Base& setHausdorff(Real hausd)
      {
        m_hausd = hausd;
        return *this;
      }

      /**
       * @brief Sets the gradation parameter
       *
       * @param[in] hgrad Gradation parameter
       *
       * @returns Reference to self (for method chaining)
       *
       * The `hgrad` option allows to set the gradation value. It controls the
       * ratio between two adjacent edges. With a gradation of @f$ h @f$, two
       * adjacent edges @f$ e_1 @f$ and @f$ e_2 @f$ must respect the following
       * constraint:
       *
       * @f[
       *  \dfrac{1}{h} \leq \dfrac{ |e_1| }{ |e_2| } \leq h
       * @f]
       *
       * By default, the gradation value is 1.3.
       *
       */
      Base& setGradation(Real hgrad)
      {
        m_hgrad = hgrad;
        return *this;
      }

    protected:
      /**
       * @brief Applies configured remeshing parameters to an MMG mesh.
       * @param[in, out] mesh Target MMG mesh.
       * @returns Reference to this object.
       */
      Base& setParameters(MMG5_pMesh mesh);

    private:
      Optional<Real> m_hmin, m_hmax, m_hausd,
        m_hgrad; ///< Optional MMG scalar parameters.
      bool m_ridgeDetection; ///< Whether ridge angle detection is enabled.
  };
}
#endif
