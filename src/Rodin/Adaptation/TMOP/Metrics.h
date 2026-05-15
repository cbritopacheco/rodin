/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TMOP_METRICS_H
#define RODIN_ADAPTATION_TMOP_METRICS_H

#include <cmath>
#include <vector>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Polytope.h"
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
   * @brief Legacy constant identity target for detached prototype code.
   */
  class IdentityTargetEvaluator
  {
    public:
      Target evaluate() const
      {
        return {};
      }
  };

  /**
   * @brief Production identity target Jacobian for the strict TMOP term.
   *
   * The strict quality term evaluates @f$T = A W^{-1}@f$ at element
   * quadrature/sample points. This target supplies @f$W = I@f$, so the
   * weighted Jacobian is the reference-coordinate Jacobian of the deformed
   * mesh coordinate map.
   */
  class IdentityTargetJacobian
  {
    public:
      Matrix2 evaluate(
          const Geometry::Polytope&,
          const Math::SpatialPoint&) const
      {
        return Matrix2::Identity();
      }
  };

  inline Matrix2 linearCellJacobian2D(const Geometry::Polytope& cell)
  {
    const auto& mesh = cell.getMesh();
    const auto& vertices = cell.getVertices();
    const auto x0 = mesh.getVertexCoordinates(vertices(0));
    const auto x1 = mesh.getVertexCoordinates(vertices(1));
    const auto x2 = mesh.getVertexCoordinates(vertices(2));

    Matrix2 W;
    W(0, 0) = x1[0] - x0[0];
    W(0, 1) = x2[0] - x0[0];
    W(1, 0) = x1[1] - x0[1];
    W(1, 1) = x2[1] - x0[1];
    return W;
  }

  /**
   * @brief Fixed per-element target Jacobian captured from an initial mesh.
   *
   * This target implements the common TMOP pattern where @f$W@f$ is frozen at
   * the start of one nonlinear solve. With this target, the strict quality
   * energy is zero at @f$u = 0@f$ if the current mesh is the same mesh used to
   * construct the target. Passing a separate reference mesh lets tests verify
   * that the residual/tangent drive a perturbed mesh back toward that fixed
   * target without changing the TMOP formulation.
   */
  class InitialElementTargetJacobian
  {
    public:
      InitialElementTargetJacobian() = default;

      template <class Mesh>
      explicit InitialElementTargetJacobian(const Mesh& mesh)
      {
        const auto& conn = mesh.getConnectivity();
        m_targets.resize(mesh.getCellCount(), Matrix2::Identity());
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          if (conn.getGeometry(2, cellIndex)
              != Geometry::Polytope::Type::Triangle)
            continue;
          auto cellIterator = mesh.getPolytope(2, cellIndex);
          m_targets[static_cast<size_t>(cellIndex)] =
            linearCellJacobian2D(*cellIterator);
        }
      }

      Matrix2 evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return linearCellJacobian2D(cell);
      }

    private:
      std::vector<Matrix2> m_targets;
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
   * Strict TMOP metric:
   *   mu(T) = 1/2 ||T - I||_F^2
   *
   * This first metric is deliberately small, but it is a real target-matrix
   * TMOP metric because it is evaluated on @f$T = A W^{-1}@f$, not on edge
   * lengths or triangle quality heuristics.
   */
  class SquaredDistanceMetric final : public MetricBase
  {
    public:
      Real value(const Matrix2& T) const override
      {
        return Real(0.5) * (T - Matrix2::Identity()).squaredNorm();
      }

      Matrix2 gradient(const Matrix2& T) const
      {
        return T - Matrix2::Identity();
      }

      Matrix2 hessianAction(const Matrix2&, const Matrix2& dT) const
      {
        return dT;
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
