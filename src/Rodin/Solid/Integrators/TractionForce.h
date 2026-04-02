/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TractionForce.h
 * @brief Traction (Neumann) boundary condition for solid mechanics.
 *
 * Represents the traction contribution to the external virtual work:
 * @f[
 *   R_{\text{ext,traction}}(\mathbf{v})
 *     = \int_{\Gamma_N} \mathbf{t} \cdot \mathbf{v} \, dS
 * @f]
 * where @f$ \mathbf{t} @f$ is the prescribed traction on the Neumann
 * boundary @f$ \Gamma_N @f$.
 *
 * @note For assembly in the Rodin framework, the user should use the standard
 * @c BoundaryIntegral(t, v).over(GammaN) syntax. This header provides
 * documentation and a thin wrapper for semantic clarity.
 */
#ifndef RODIN_SOLID_INTEGRATORS_TRACTIONFORCE_H
#define RODIN_SOLID_INTEGRATORS_TRACTIONFORCE_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Solid
{
  /**
   * @brief Describes a traction (surface force) boundary condition.
   *
   * In a Rodin variational formulation:
   * @code
   * VectorFunction t{1.0, 0.0};
   * problem = ...
   *         - BoundaryIntegral(t, v).over(GammaN);
   * @endcode
   *
   * The TractionForce class provides semantic clarity:
   * @code
   * TractionForce traction({1.0, 0.0});
   * // traction.getValue() returns the traction vector
   * @endcode
   */
  class TractionForce
  {
    public:
      /**
       * @brief Constructs a traction force from a constant vector.
       * @param traction Traction vector @f$ \mathbf{t} @f$
       */
      TractionForce(const Math::SpatialVector<Real>& traction)
        : m_traction(traction)
      {}

      TractionForce(const TractionForce&) = default;
      TractionForce(TractionForce&&) = default;

      /// @brief Gets the traction vector.
      const Math::SpatialVector<Real>& getValue() const { return m_traction; }

      /// @brief Sets the traction vector.
      TractionForce& setValue(const Math::SpatialVector<Real>& traction)
      {
        m_traction = traction;
        return *this;
      }

    private:
      Math::SpatialVector<Real> m_traction;
  };
}

#endif
