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
        m_H.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        m_H.setZero();
        m_HDense.resize(d, d);
        m_HDense.setZero();
        m_F.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        m_F.setIdentity();
        m_FDense.resize(d, d);
        m_FDense.setIdentity();
        m_Finv.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        m_Finv.setIdentity();
        m_FinvDense.resize(d, d);
        m_FinvDense.setIdentity();
        m_FinvT.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        m_FinvT.setIdentity();
        m_FinvTDense.resize(d, d);
        m_FinvTDense.setIdentity();
        m_C.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        m_C.setIdentity();
        m_CDense.resize(d, d);
        m_CDense.setIdentity();
        m_b.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        m_b.setIdentity();
        m_bDense.resize(d, d);
        m_bDense.setIdentity();
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
        Math::SpatialMatrix<Real> HSpatial;
        HSpatial = H;
        return setDisplacementGradient(HSpatial);
      }

      /**
       * @brief Sets the displacement gradient @f$ \nabla \mathbf{u} @f$.
       * @param H Displacement gradient matrix
       * @returns Reference to this for chaining
       */
      KinematicState& setDisplacementGradient(const Math::SpatialMatrix<Real>& H)
      {
        updateKinematics(H);
        return *this;
      }

      /// @brief Gets the spatial dimension.
      size_t getDimension() const { return m_d; }

      /// @brief Gets the displacement gradient @f$ \nabla \mathbf{u} @f$.
      const Math::Matrix<Real>& getDisplacementGradient() const { return m_HDense; }

      /// @brief Gets the deformation gradient @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$.
      const Math::Matrix<Real>& getDeformationGradient() const { return m_FDense; }

      /// @brief Gets @f$ \mathbf{F}^{-1} @f$.
      const Math::Matrix<Real>& getDeformationGradientInverse() const { return m_FinvDense; }

      /// @brief Gets @f$ \mathbf{F}^{-T} @f$.
      const Math::Matrix<Real>& getDeformationGradientInverseTranspose() const { return m_FinvTDense; }

      /// @brief Gets the right Cauchy-Green tensor @f$ \mathbf{C} = \mathbf{F}^T \mathbf{F} @f$.
      const Math::Matrix<Real>& getRightCauchyGreenTensor() const { return m_CDense; }

      /// @brief Gets the left Cauchy-Green tensor @f$ \mathbf{b} = \mathbf{F} \mathbf{F}^T @f$.
      const Math::Matrix<Real>& getLeftCauchyGreenTensor() const { return m_bDense; }

      /// @brief Gets the Jacobian @f$ J = \det(\mathbf{F}) @f$.
      Real getJacobian() const { return m_J; }

      /// @brief Gets @f$ \ln(J) @f$.
      Real getLogJacobian() const { return m_logJ; }

    private:
      size_t m_d;                  ///< Spatial dimension

      Math::SpatialMatrix<Real> m_H;     ///< Displacement gradient
      Math::SpatialMatrix<Real> m_F;     ///< Deformation gradient
      Math::SpatialMatrix<Real> m_Finv;  ///< Inverse of deformation gradient
      Math::SpatialMatrix<Real> m_FinvT; ///< Inverse transpose of deformation gradient
      Math::SpatialMatrix<Real> m_C;     ///< Right Cauchy-Green tensor
      Math::SpatialMatrix<Real> m_b;     ///< Left Cauchy-Green tensor

      Math::Matrix<Real> m_HDense;     ///< Dense displacement gradient (API compatibility)
      Math::Matrix<Real> m_FDense;     ///< Dense deformation gradient (API compatibility)
      Math::Matrix<Real> m_FinvDense;  ///< Dense inverse deformation gradient (API compatibility)
      Math::Matrix<Real> m_FinvTDense; ///< Dense inverse transpose deformation gradient (API compatibility)
      Math::Matrix<Real> m_CDense;     ///< Dense right Cauchy-Green tensor (API compatibility)
      Math::Matrix<Real> m_bDense;     ///< Dense left Cauchy-Green tensor (API compatibility)
      Real m_J;                    ///< Jacobian (det F)
      Real m_logJ;                 ///< log(J)

      static
      void spatialToDense(Math::Matrix<Real>& dst, const Math::SpatialMatrix<Real>& src)
      {
        // Extract the active rows/cols from bounded SpatialMatrix storage into
        // a dynamic dense matrix for API compatibility.
        const size_t rows = static_cast<size_t>(src.rows());
        const size_t cols = static_cast<size_t>(src.cols());
        dst.resize(rows, cols);
        dst = src.getData().topLeftCorner(
          static_cast<Eigen::Index>(rows),
          static_cast<Eigen::Index>(cols));
      }

      void syncDenseViews()
      {
        // Keep dense views synchronized with internal spatial storage so
        // existing callers using Math::Matrix<Real> getters remain unchanged.
        spatialToDense(m_HDense, m_H);
        spatialToDense(m_FDense, m_F);
        spatialToDense(m_FinvDense, m_Finv);
        spatialToDense(m_FinvTDense, m_FinvT);
        spatialToDense(m_CDense, m_C);
        spatialToDense(m_bDense, m_b);
      }

      void updateKinematics(const Math::SpatialMatrix<Real>& H)
      {
        assert(static_cast<size_t>(H.rows()) == m_d);
        assert(static_cast<size_t>(H.cols()) == m_d);
        m_H = H;

        // F = I + H
        m_F = m_H;
        for (size_t i = 0; i < m_d; i++)
          m_F(static_cast<std::uint8_t>(i), static_cast<std::uint8_t>(i)) += 1.0;

        // J = det(F)
        m_J = m_F.determinant();

        // F^{-1} and F^{-T}
        m_Finv = m_F.inverse();
        m_FinvT = m_Finv.transpose();

        // C = F^T F
        m_C = m_F.transpose() * m_F;

        // b = F F^T
        m_b = m_F * m_F.transpose();

        syncDenseViews();

        // log(J)
        m_logJ = std::log(m_J);
      }
  };
}

#endif
