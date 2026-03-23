/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Invariants.h
 * @brief Strain invariants for isotropic and anisotropic hyperelasticity.
 *
 * Provides isotropic invariants @f$ I_1, I_2, I_3 @f$ of the right
 * Cauchy-Green tensor @f$ \mathbf{C} @f$:
 * @f[
 *   I_1 = \operatorname{tr}(\mathbf{C}), \quad
 *   I_2 = \tfrac{1}{2}[\operatorname{tr}(\mathbf{C})^2 - \operatorname{tr}(\mathbf{C}^2)], \quad
 *   I_3 = \det(\mathbf{C}) = J^2
 * @f]
 *
 * and fiber invariants @f$ I_4, I_5 @f$ for structural tensors:
 * @f[
 *   I_4 = \mathbf{a}_0 \cdot \mathbf{C} \, \mathbf{a}_0, \quad
 *   I_5 = \mathbf{a}_0 \cdot \mathbf{C}^2 \, \mathbf{a}_0
 * @f]
 */
#ifndef RODIN_SOLID_KINEMATICS_INVARIANTS_H
#define RODIN_SOLID_KINEMATICS_INVARIANTS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"

#include "KinematicState.h"

namespace Rodin::Solid
{
  /**
   * @brief Isotropic invariants of the right Cauchy-Green tensor.
   *
   * Computes the three principal invariants:
   * - @f$ I_1 = \operatorname{tr}(\mathbf{C}) @f$
   * - @f$ I_2 = \tfrac{1}{2}[\operatorname{tr}(\mathbf{C})^2 - \operatorname{tr}(\mathbf{C}^2)] @f$
   * - @f$ I_3 = \det(\mathbf{C}) = J^2 @f$
   *
   * ## Usage
   * @code
   * IsotropicInvariants inv;
   * inv.setState(state);
   * Real I1 = inv.getFirstInvariant();
   * @endcode
   */
  class IsotropicInvariants
  {
    public:
      IsotropicInvariants() = default;

      /**
       * @brief Sets the kinematic state and computes all invariants.
       * @param state Kinematic state from which @f$ \mathbf{C} @f$ and @f$ J @f$ are read
       * @returns Reference to this for chaining
       */
      IsotropicInvariants& setState(const KinematicState& state)
      {
        const auto& C = state.getRightCauchyGreenTensor();
        const Real J = state.getJacobian();

        m_I1 = C.trace();
        m_I2 = 0.5 * (m_I1 * m_I1 - (C * C).trace());
        m_I3 = J * J;

        return *this;
      }

      /// @brief Gets @f$ I_1 = \operatorname{tr}(\mathbf{C}) @f$.
      Real getFirstInvariant() const { return m_I1; }

      /// @brief Gets @f$ I_2 = \tfrac{1}{2}[\operatorname{tr}(\mathbf{C})^2 - \operatorname{tr}(\mathbf{C}^2)] @f$.
      Real getSecondInvariant() const { return m_I2; }

      /// @brief Gets @f$ I_3 = \det(\mathbf{C}) = J^2 @f$.
      Real getThirdInvariant() const { return m_I3; }

    private:
      Real m_I1 = 0.0;
      Real m_I2 = 0.0;
      Real m_I3 = 0.0;
  };

  /**
   * @brief Fiber invariants for anisotropic hyperelasticity.
   *
   * Computes invariants associated with a preferred fiber direction
   * @f$ \mathbf{a}_0 @f$:
   * - @f$ I_4 = \mathbf{a}_0 \cdot \mathbf{C} \, \mathbf{a}_0 @f$
   * - @f$ I_5 = \mathbf{a}_0 \cdot \mathbf{C}^2 \, \mathbf{a}_0 @f$
   *
   * ## Usage
   * @code
   * Math::SpatialVector<Real> fiberDir(2);
   * fiberDir << 1.0, 0.0;
   * FiberInvariants inv(fiberDir);
   * inv.setState(state);
   * Real I4 = inv.getFourthInvariant();
   * @endcode
   */
  class FiberInvariants
  {
    public:
      /**
       * @brief Constructs fiber invariants for a given fiber direction.
       * @param a0 Unit fiber direction vector @f$ \mathbf{a}_0 @f$
       */
      FiberInvariants(const Math::SpatialVector<Real>& a0)
        : m_a0(a0)
      {}

      FiberInvariants(const FiberInvariants&) = default;
      FiberInvariants(FiberInvariants&&) = default;

      /**
       * @brief Sets the kinematic state and computes fiber invariants.
       * @param state Kinematic state from which @f$ \mathbf{C} @f$ is read
       * @returns Reference to this for chaining
       */
      FiberInvariants& setState(const KinematicState& state)
      {
        const auto& C = state.getRightCauchyGreenTensor();
        m_I4 = m_a0.dot(C * m_a0);
        m_I5 = m_a0.dot(C * C * m_a0);
        return *this;
      }

      /// @brief Gets @f$ I_4 = \mathbf{a}_0 \cdot \mathbf{C} \, \mathbf{a}_0 @f$.
      Real getFourthInvariant() const { return m_I4; }

      /// @brief Gets @f$ I_5 = \mathbf{a}_0 \cdot \mathbf{C}^2 \, \mathbf{a}_0 @f$.
      Real getFifthInvariant() const { return m_I5; }

    private:
      Math::SpatialVector<Real> m_a0;
      Real m_I4 = 0.0;
      Real m_I5 = 0.0;
  };
}

#endif
