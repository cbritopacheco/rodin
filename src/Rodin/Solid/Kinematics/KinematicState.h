/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file KinematicState.h
 * @brief Kinematic state for finite-strain solid mechanics.
 *
 * Provides the KinematicState class that encapsulates all kinematic quantities
 * derived from the displacement gradient, including:
 * - Deformation gradient @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$
 * - Right Cauchy-Green tensor @f$ \mathbf{C} = \mathbf{F}^T \mathbf{F} @f$
 * - Left Cauchy-Green tensor @f$ \mathbf{b} = \mathbf{F} \mathbf{F}^T @f$
 * - Jacobian @f$ J = \det(\mathbf{F}) @f$
 * - Inverse and inverse transpose of @f$ \mathbf{F} @f$
 */
#ifndef RODIN_SOLID_KINEMATICS_KINEMATICSTATE_H
#define RODIN_SOLID_KINEMATICS_KINEMATICSTATE_H

#include <cassert>
#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Solid
{
  /**
   * @brief Kinematic state for finite-strain continuum mechanics.
   *
   * Stores the displacement gradient and computes all standard kinematic
   * quantities needed for hyperelastic formulations:
   *
   * | Quantity | Symbol | Definition |
   * |----------|--------|------------|
   * | Deformation gradient | @f$ \mathbf{F} @f$ | @f$ \mathbf{I} + \nabla \mathbf{u} @f$ |
   * | Right Cauchy-Green | @f$ \mathbf{C} @f$ | @f$ \mathbf{F}^T \mathbf{F} @f$ |
   * | Left Cauchy-Green | @f$ \mathbf{b} @f$ | @f$ \mathbf{F} \mathbf{F}^T @f$ |
   * | Jacobian | @f$ J @f$ | @f$ \det(\mathbf{F}) @f$ |
   *
   * ## Usage
   * @code
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU);
   * state.update();
   * auto F = state.getDeformationGradient();
   * Real J = state.getJacobian();
   * @endcode
   */
  class KinematicState
  {
    public:
      /**
       * @brief Constructs a kinematic state for the given spatial dimension.
       * @param d Spatial dimension (2 or 3)
       */
      KinematicState(size_t d)
        : m_d(d)
      {
        assert(d == 2 || d == 3);
        m_H.resize(d, d);
        m_H.setZero();
        m_F.resize(d, d);
        m_F.setIdentity();
        m_Finv.resize(d, d);
        m_Finv.setIdentity();
        m_FinvT.resize(d, d);
        m_FinvT.setIdentity();
        m_C.resize(d, d);
        m_C.setIdentity();
        m_b.resize(d, d);
        m_b.setIdentity();
        m_J = 1.0;
        m_logJ = 0.0;
      }

      KinematicState(const KinematicState&) = default;
      KinematicState(KinematicState&&) = default;
      KinematicState& operator=(const KinematicState&) = default;
      KinematicState& operator=(KinematicState&&) = default;

      /**
       * @brief Sets the displacement gradient @f$ \nabla \mathbf{u} @f$.
       * @param H Displacement gradient matrix
       * @returns Reference to this for chaining
       */
      KinematicState& setDisplacementGradient(const Math::Matrix<Real>& H)
      {
        assert(static_cast<size_t>(H.rows()) == m_d);
        assert(static_cast<size_t>(H.cols()) == m_d);
        m_H = H;
        return *this;
      }

      /**
       * @brief Updates all derived kinematic quantities from the current
       * displacement gradient.
       *
       * Computes @f$ \mathbf{F}, \mathbf{C}, \mathbf{b}, J @f$ and
       * the inverse/transpose quantities.
       *
       * @returns Reference to this for chaining
       */
      KinematicState& update()
      {
        // F = I + H
        m_F = Math::Matrix<Real>::Identity(m_d, m_d) + m_H;

        // J = det(F)
        m_J = m_F.determinant();

        // F^{-1} and F^{-T}
        m_Finv = m_F.inverse();
        m_FinvT = m_Finv.transpose();

        // C = F^T F
        m_C = m_F.transpose() * m_F;

        // b = F F^T
        m_b = m_F * m_F.transpose();

        // log(J)
        m_logJ = std::log(m_J);

        return *this;
      }

      /// @brief Gets the spatial dimension.
      size_t getDimension() const { return m_d; }

      /// @brief Gets the displacement gradient @f$ \nabla \mathbf{u} @f$.
      const Math::Matrix<Real>& getDisplacementGradient() const { return m_H; }

      /// @brief Gets the deformation gradient @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$.
      const Math::Matrix<Real>& getDeformationGradient() const { return m_F; }

      /// @brief Gets @f$ \mathbf{F}^{-1} @f$.
      const Math::Matrix<Real>& getDeformationGradientInverse() const { return m_Finv; }

      /// @brief Gets @f$ \mathbf{F}^{-T} @f$.
      const Math::Matrix<Real>& getDeformationGradientInverseTranspose() const { return m_FinvT; }

      /// @brief Gets the right Cauchy-Green tensor @f$ \mathbf{C} = \mathbf{F}^T \mathbf{F} @f$.
      const Math::Matrix<Real>& getRightCauchyGreenTensor() const { return m_C; }

      /// @brief Gets the left Cauchy-Green tensor @f$ \mathbf{b} = \mathbf{F} \mathbf{F}^T @f$.
      const Math::Matrix<Real>& getLeftCauchyGreenTensor() const { return m_b; }

      /// @brief Gets the Jacobian @f$ J = \det(\mathbf{F}) @f$.
      Real getJacobian() const { return m_J; }

      /// @brief Gets @f$ \ln(J) @f$.
      Real getLogJacobian() const { return m_logJ; }

    private:
      size_t m_d;                  ///< Spatial dimension

      Math::Matrix<Real> m_H;     ///< Displacement gradient
      Math::Matrix<Real> m_F;     ///< Deformation gradient
      Math::Matrix<Real> m_Finv;  ///< Inverse of deformation gradient
      Math::Matrix<Real> m_FinvT; ///< Inverse transpose of deformation gradient
      Math::Matrix<Real> m_C;     ///< Right Cauchy-Green tensor
      Math::Matrix<Real> m_b;     ///< Left Cauchy-Green tensor
      Real m_J;                    ///< Jacobian (det F)
      Real m_logJ;                 ///< log(J)
  };
}

#endif
