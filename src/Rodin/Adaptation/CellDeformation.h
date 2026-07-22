/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_CELLDEFORMATION_H
#define RODIN_ADAPTATION_CELLDEFORMATION_H

#include <cmath>

#include "Rodin/Math.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Lazily evaluated deformation state of a mesh cell.
   *
   * Given the displacement gradient @f$H=\nabla u@f$ at a point, this holds the
   * deformation gradient @f$F=I+H@f$ and derives, on demand, the quantities mesh
   * adaptation needs: the Jacobian @f$j=\det F@f$, the inverse transpose
   * @f$F^{-\top}@f$, the relative distortion
   * @f[
   *   Q_{\operatorname{rel}}(F)=\frac{\lVert F\rVert_F^2}{d\,j^{2/d}},
   * @f]
   * and its gradient @f$\partial Q_{\operatorname{rel}}/\partial F@f$. Each is
   * computed at most once per state and cached, following the lazy
   * @c mutable @ref Rodin::Optional idiom of @ref Rodin::Geometry::Point.
   *
   * This is the mesh-adaptation counterpart of
   * @ref Rodin::Solid::Kinematics::KinematicState. It differs deliberately in two
   * ways, both required here:
   * - it is *lazy*: an admissibility sweep needs only @f$j@f$ at most points, and
   *   never needs @f$C@f$, @f$b@f$ or @f$\log j@f$, so computing them eagerly
   *   would be wasted work at every validation point;
   * - it does *not* assume @f$j>0@f$. Mesh adaptation must evaluate inverted and
   *   degenerate states precisely in order to detect and reject them, so
   *   admissibility is reported through @ref isAdmissible rather than asserted.
   *   Quantities that are undefined for @f$j\le0@f$ (@ref getInverseTranspose,
   *   @ref getRelativeDistortion, @ref getRelativeDistortionGradient) must only be
   *   requested once @ref isAdmissible holds.
   */
  class CellDeformation
  {
    public:
      /// @brief Constructs an undeformed state of the given spatial dimension.
      explicit CellDeformation(std::size_t d)
        : m_d(d)
      {
        assert(d == 2 || d == 3);
        m_F = Math::SpatialMatrix<Real>::Identity(
          static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
      }

      /**
       * @brief Sets the displacement gradient @f$H=\nabla u@f$, giving
       * @f$F=I+H@f$, and invalidates the derived quantities.
       */
      CellDeformation& setDisplacementGradient(
        const Math::SpatialMatrix<Real>& H)
      {
        m_F = Math::SpatialMatrix<Real>::Identity(
          static_cast<std::uint8_t>(m_d), static_cast<std::uint8_t>(m_d));
        m_F += H;
        invalidate();
        return *this;
      }

      /// @brief Sets the deformation gradient @f$F@f$ directly.
      CellDeformation& setDeformationGradient(
        const Math::SpatialMatrix<Real>& F)
      {
        m_F = F;
        invalidate();
        return *this;
      }

      /// @brief The spatial dimension.
      std::size_t getDimension() const { return m_d; }

      /// @brief The deformation gradient @f$F=I+\nabla u@f$.
      const Math::SpatialMatrix<Real>& getDeformationGradient() const
      {
        return m_F;
      }

      /// @brief The Jacobian @f$j=\det F@f$, computed once and cached.
      Real getJacobian() const
      {
        if (!m_j)
          m_j = m_F.determinant();
        return *m_j;
      }

      /**
       * @brief Whether the cell is non-inverted, @f$j>0@f$.
       *
       * The derived shape quantities are defined only in this case.
       */
      bool isAdmissible() const { return getJacobian() > Real(0); }

      /**
       * @brief Whether @f$F@f$ is numerically invertible, @f$|j|>\varepsilon@f$.
       *
       * Weaker than @ref isAdmissible: an inverted cell (@f$j<0@f$) is
       * invertible. Mesh adaptation relies on this distinction, since the
       * size-control terms act on inverted cells precisely to pull them back to
       * validity, and need @f$F^{-\top}@f$ there.
       */
      bool isInvertible() const
      {
        return std::abs(getJacobian()) > s_singularTolerance;
      }

      /**
       * @brief The inverse transpose @f$F^{-\top}@f$; requires
       * @ref isInvertible.
       *
       * Defined for inverted cells too: only the cofactor structure is needed,
       * not the sign of @f$j@f$.
       */
      const Math::SpatialMatrix<Real>& getInverseTranspose() const
      {
        assert(isInvertible());
        if (!m_FinvT)
          m_FinvT = Math::SpatialMatrix<Real>(m_F.inverse().transpose());
        return *m_FinvT;
      }

      /**
       * @brief The relative distortion @f$Q_{\operatorname{rel}}@f$; requires
       * @ref isAdmissible.
       *
       * Invariant under @f$F\mapsto sRF@f$ for scalars @f$s>0@f$ and rotations
       * @f$R@f$, and minimal at similarities.
       */
      Real getRelativeDistortion() const
      {
        assert(isAdmissible());
        if (!m_Q)
        {
          const Real d = static_cast<Real>(m_d);
          m_Q = m_F.squaredNorm() /
            (d * std::pow(getJacobian(), Real(2) / d));
        }
        return *m_Q;
      }

      /**
       * @brief The gradient @f$\partial Q_{\operatorname{rel}}/\partial F@f$;
       * requires @ref isAdmissible.
       */
      const Math::SpatialMatrix<Real>& getRelativeDistortionGradient() const
      {
        assert(isAdmissible());
        if (!m_dQdF)
        {
          const Real d = static_cast<Real>(m_d);
          const Real j = getJacobian();
          const Real frob2 = m_F.squaredNorm();
          m_dQdF = Math::SpatialMatrix<Real>(
            (Real(2) / d) * std::pow(j, -Real(2) / d) *
            (m_F - (frob2 / d) * getInverseTranspose()));
        }
        return *m_dQdF;
      }

      /**
       * @brief The linearised action of the Jacobian on a gradient perturbation,
       * @f$a_j(G)=\partial_F(\det F):G=j\,F^{-\top}\!:G@f$.
       *
       * This is the directional derivative of @f$j@f$ along @f$G=\nabla v@f$,
       * and requires only @ref isInvertible.
       */
      Real getJacobianAction(const Math::SpatialMatrix<Real>& G) const
      {
        return getJacobian() * getInverseTranspose().dot(G);
      }

      /**
       * @brief The linearised action of the relative distortion,
       * @f$a_Q(G)=\partial_F Q_{\operatorname{rel}}:G@f$; requires
       * @ref isAdmissible.
       */
      Real getRelativeDistortionAction(const Math::SpatialMatrix<Real>& G) const
      {
        return getRelativeDistortionGradient().dot(G);
      }

    private:
      /// @brief Below this @f$|j|@f$ the cofactor via the inverse is unreliable.
      static constexpr Real s_singularTolerance = Real(1e-14);

      void invalidate() const
      {
        m_j.reset();
        m_FinvT.reset();
        m_Q.reset();
        m_dQdF.reset();
      }

      std::size_t m_d;
      Math::SpatialMatrix<Real> m_F;
      mutable Optional<Real> m_j;
      mutable Optional<Real> m_Q;
      mutable Optional<Math::SpatialMatrix<Real>> m_FinvT;
      mutable Optional<Math::SpatialMatrix<Real>> m_dQdF;
  };
}

#endif
