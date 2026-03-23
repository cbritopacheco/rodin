/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BodyForce.h
 * @brief Body force contribution for solid mechanics formulations.
 *
 * Represents the body force contribution to the external virtual work:
 * @f[
 *   R_{\text{ext,body}}(\mathbf{v}) = \int_\Omega \mathbf{f} \cdot \mathbf{v} \, dX
 * @f]
 * where @f$ \mathbf{f} @f$ is the body force per unit reference volume.
 *
 * @note For assembly in the Rodin framework, the user should use the standard
 * @c Integral(f, v) syntax. This header provides documentation and a
 * thin wrapper for semantic clarity in a solid mechanics context.
 */
#ifndef RODIN_SOLID_INTEGRATORS_BODYFORCE_H
#define RODIN_SOLID_INTEGRATORS_BODYFORCE_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Solid
{
  /**
   * @brief Describes a body force per unit reference volume.
   *
   * This is a lightweight container for body force data. In a Rodin
   * variational formulation, the user assembles the load vector using:
   * @code
   * VectorFunction f{0.0, -9.81};
   * problem = ...
   *         - Integral(f, v);
   * @endcode
   *
   * The BodyForce class serves as a semantic wrapper when working
   * within the Solid module's API:
   * @code
   * BodyForce body({0.0, -9.81});
   * // body.getValue() returns the force vector
   * @endcode
   */
  class BodyForce
  {
    public:
      /**
       * @brief Constructs a body force from a constant vector.
       * @param force Body force vector @f$ \mathbf{f} @f$
       */
      BodyForce(const Math::SpatialVector<Real>& force)
        : m_force(force)
      {}

      BodyForce(const BodyForce&) = default;
      BodyForce(BodyForce&&) = default;

      /// @brief Gets the body force vector.
      const Math::SpatialVector<Real>& getValue() const { return m_force; }

      /// @brief Sets the body force vector.
      BodyForce& setValue(const Math::SpatialVector<Real>& force)
      {
        m_force = force;
        return *this;
      }

    private:
      Math::SpatialVector<Real> m_force;
  };
}

#endif
