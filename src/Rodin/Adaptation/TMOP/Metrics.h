/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TMOP_METRICS_H
#define RODIN_ADAPTATION_TMOP_METRICS_H

#include <cmath>

#include "Rodin/Math.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation::TMOP
{
  using Matrix2 = Math::FixedSizeMatrix<Real, 2, 2>;

  struct Target
  {
    Matrix2 W = Matrix2::Identity();
  };

  /**
   * @brief Constant identity target for early TMOP infrastructure tests.
   *
   * Future target evaluators can supply element- and quadrature-dependent
   * targets. This first one keeps the weighted Jacobian equal to the geometric
   * Jacobian, which is sufficient for exercising the metric/objective pipeline.
   */
  class IdentityTargetEvaluator
  {
    public:
      Target evaluate() const
      {
        return {};
      }
  };

  class MetricBase
  {
    public:
      virtual ~MetricBase() = default;

      virtual Real value(const Matrix2& A) const = 0;

      /**
       * @brief Sets the determinant threshold below which a Jacobian is invalid.
       */
      MetricBase& setDeterminantTolerance(Real eps)
      {
        m_detTolerance = std::max(Real(0), eps);
        return *this;
      }

      MetricBase& setInvalidPenalty(Real penalty)
      {
        m_invalidPenalty = penalty;
        return *this;
      }

      Real getDeterminantTolerance() const
      {
        return m_detTolerance;
      }

      Real getInvalidPenalty() const
      {
        return m_invalidPenalty;
      }

    protected:
      /**
       * @brief Returns zero for admissible Jacobians and a large penalty otherwise.
       *
       * The penalty is intentionally simple. It guards objective evaluation
       * against inverted sampled Jacobians.
       */
      Real invalidPenalty(const Matrix2& A) const
      {
        const Real det = A.determinant();
        if (det > m_detTolerance)
          return Real(0);
        return m_invalidPenalty + (m_detTolerance - det) * m_invalidPenalty;
      }

    private:
      Real m_detTolerance = 1e-12;
      Real m_invalidPenalty = 1e20;
  };

  /**
   * Simple infrastructure metric:
   *   W(A) = 1/2 ||A - I||_F^2
   *
   * This is the initial differentiable quality metric used by the TMOP
   * evaluator and residual/tangent terms.
   */
  class SquaredDistanceMetric final : public MetricBase
  {
    public:
      Real value(const Matrix2& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        return Real(0.5) * (A - Matrix2::Identity()).squaredNorm();
      }
  };

  /**
   * Frobenius shape distortion in 2D:
   *   W(A) = ||A||_F^2 / (2 det(A)) - 1
   *
   * The metric is zero on rotations/scaled-orthogonal maps with unit shape and
   * grows near singularity. Non-positive determinants receive a large penalty.
   */
  class ShapeDistortionMetric final : public MetricBase
  {
    public:
      Real value(const Matrix2& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        return A.squaredNorm() / (Real(2) * A.determinant()) - Real(1);
      }
  };

  /**
   * Area distortion metric:
   *   W(A) = 1/2 (det(A) - 1)^2
   */
  class AreaDistortionMetric final : public MetricBase
  {
    public:
      Real value(const Matrix2& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        const Real d = A.determinant() - Real(1);
        return Real(0.5) * d * d;
      }
  };

  class ShapeSizeMetric final : public MetricBase
  {
    public:
      Real value(const Matrix2& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        return m_shape.value(A) + m_area.value(A);
      }

    private:
      ShapeDistortionMetric m_shape;
      AreaDistortionMetric m_area;
  };
}

#endif
