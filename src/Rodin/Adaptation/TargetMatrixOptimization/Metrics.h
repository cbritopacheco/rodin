/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_METRICS_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_METRICS_H

#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
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

  /**
   * @brief Size-consistent equilateral target Jacobian.
   *
   * For each triangle @f$K@f$ this captures the initial element area
   * @f$|K_0|@f$ and sets @f$W_K@f$ to the reference-to-physical Jacobian of an
   * equilateral triangle with that same area:
   * @f[
   *   l_K = \sqrt{4|K_0|/\sqrt 3},\qquad
   *   W_K = \begin{pmatrix} l_K & l_K/2 \\ 0 & l_K\sqrt3/2 \end{pmatrix},
   *   \qquad \tfrac12\det(W_K) = |K_0|.
   * @f]
   *
   * Because @f$\det(W_K)@f$ is in the physical element scale, @f$T = A W^{-1}@f$
   * is @f$O(1)@f$ even on a domain scaled to @f$[0,1]@f$, so the quality term
   * recovers element *shape* toward equilateral without trying to rescale every
   * element to unit size. This is the correct production target;
   * IdentityTargetJacobian (@f$W=I@f$) is only a unit-scale sanity target.
   */
  class EquilateralSameAreaTargetJacobian
  {
    public:
      EquilateralSameAreaTargetJacobian() = default;

      template <class Mesh>
      explicit EquilateralSameAreaTargetJacobian(const Mesh& mesh)
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
            equilateralSameArea(linearCellJacobian2D(*cellIterator));
        }
      }

      Matrix2 evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return equilateralSameArea(linearCellJacobian2D(cell));
      }

      static Matrix2 equilateralSameArea(const Matrix2& A0)
      {
        const Real area = Real(0.5) * std::abs(A0.determinant());
        // Degenerate initial element: keep its (tiny) Jacobian so W stays in
        // physical scale instead of collapsing to a singular target.
        if (area <= Real(1e-30))
          return A0;
        const Real l = std::sqrt(Real(4) * area / std::sqrt(Real(3)));
        Matrix2 W;
        W(0, 0) = l;
        W(0, 1) = Real(0.5) * l;
        W(1, 0) = Real(0);
        W(1, 1) = l * std::sqrt(Real(3)) / Real(2);
        return W;
      }

    private:
      std::vector<Matrix2> m_targets;
  };

  /**
   * @brief Generic, polytope- and FES-independent reference Jacobian helper.
   *
   * Uses the cell's isoparametric map evaluated at the reference-element
   * centroid, so it works for any 2D polytope (triangle, quadrilateral, ...)
   * and any element order without baking in P1/triangle vertex indexing.
   */
  inline Matrix2 isoparametricCellJacobian2D(const Geometry::Polytope& cell)
  {
    const Geometry::Polytope::Traits traits(cell.getGeometry());
    const Math::SpatialPoint rc = traits.getCentroid();
    const Geometry::Point point(cell, rc);
    const auto& J = point.getJacobian();
    Matrix2 W;
    W(0, 0) = J(0, 0);
    W(0, 1) = J(0, 1);
    W(1, 0) = J(1, 0);
    W(1, 1) = J(1, 1);
    return W;
  }

  /**
   * @brief Constant target Jacobian: @f$W@f$ is a fixed matrix, independent of
   * the cell, the polytope type, and the finite element space.
   *
   * This is the most general analytic target. The factory helpers cover the
   * common TMOP cases: a uniform resize, an axis-aligned anisotropic target, a
   * pure-rotation (shape-only) target, and an oriented anisotropic stretch
   * @f$W = R(\theta)\,\mathrm{diag}(s_x, s_y)@f$.
   */
  class ConstantTargetJacobian
  {
    public:
      ConstantTargetJacobian()
        : m_W(Matrix2::Identity())
      {}

      explicit ConstantTargetJacobian(const Matrix2& W)
        : m_W(W)
      {}

      Matrix2 evaluate(
          const Geometry::Polytope&,
          const Math::SpatialPoint&) const
      {
        return m_W;
      }

      const Matrix2& getMatrix() const
      {
        return m_W;
      }

      static ConstantTargetJacobian identity()
      {
        return ConstantTargetJacobian(Matrix2::Identity());
      }

      static ConstantTargetJacobian uniformScale(Real s)
      {
        Matrix2 W = Matrix2::Identity();
        W(0, 0) = s;
        W(1, 1) = s;
        return ConstantTargetJacobian(W);
      }

      static ConstantTargetJacobian diagonal(Real sx, Real sy)
      {
        Matrix2 W = Matrix2::Zero();
        W(0, 0) = sx;
        W(1, 1) = sy;
        return ConstantTargetJacobian(W);
      }

      static ConstantTargetJacobian rotation(Real theta)
      {
        const Real c = std::cos(theta);
        const Real s = std::sin(theta);
        Matrix2 W;
        W(0, 0) = c;
        W(0, 1) = -s;
        W(1, 0) = s;
        W(1, 1) = c;
        return ConstantTargetJacobian(W);
      }

      static ConstantTargetJacobian stretch(Real sx, Real sy, Real theta)
      {
        const Real c = std::cos(theta);
        const Real s = std::sin(theta);
        // W = R(theta) * diag(sx, sy).
        Matrix2 W;
        W(0, 0) = c * sx;
        W(0, 1) = -s * sy;
        W(1, 0) = s * sx;
        W(1, 1) = c * sy;
        return ConstantTargetJacobian(W);
      }

    private:
      Matrix2 m_W;
  };

  /**
   * @brief Spatially varying target driven by a user callable
   * @f$W = f(x)@f$ of the physical point.
   *
   * Polytope- and FES-independent: it only maps the reference sample point to
   * physical coordinates through the cell's isoparametric map and evaluates the
   * callable, so it works for any element type. Use makeAnalyticTargetJacobian
   * for class template argument deduction from a lambda.
   */
  template <class Function>
  class AnalyticTargetJacobian
  {
    public:
      explicit AnalyticTargetJacobian(Function f)
        : m_f(std::move(f))
      {}

      Matrix2 evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint& rc) const
      {
        const Geometry::Point point(cell, rc);
        return m_f(point.getPhysicalCoordinates());
      }

    private:
      Function m_f;
  };

  template <class Function>
  AnalyticTargetJacobian<std::decay_t<Function>>
  makeAnalyticTargetJacobian(Function&& f)
  {
    return AnalyticTargetJacobian<std::decay_t<Function>>(
        std::forward<Function>(f));
  }

  /**
   * @brief Polytope- and FES-independent generalization of
   * InitialElementTargetJacobian.
   *
   * Captures, per cell, the initial isoparametric Jacobian evaluated at the
   * reference centroid via Geometry::Point, so it is correct for triangles,
   * quadrilaterals and higher-order elements alike. On affine triangles it
   * reproduces InitialElementTargetJacobian exactly (and thus gives zero strict
   * energy at @f$u = 0@f$ on the captured mesh), but it does not assume P1
   * triangle vertex indexing.
   */
  class IsoparametricTargetJacobian
  {
    public:
      IsoparametricTargetJacobian() = default;

      template <class Mesh>
      explicit IsoparametricTargetJacobian(const Mesh& mesh)
      {
        const auto& conn = mesh.getConnectivity();
        m_targets.resize(mesh.getCellCount(), Matrix2::Identity());
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          if (conn.getGeometry(2, cellIndex) == Geometry::Polytope::Type::Point)
            continue;
          auto cellIterator = mesh.getPolytope(2, cellIndex);
          m_targets[static_cast<size_t>(cellIndex)] =
            isoparametricCellJacobian2D(*cellIterator);
        }
      }

      Matrix2 evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return isoparametricCellJacobian2D(cell);
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
   * Frobenius shape distortion in 2D (MFEM TMOP metric 2):
   *   mu(T) = ||T||_F^2 / (2 det(T)) - 1
   *
   * Scale-invariant barrier shape metric: mu(sT) = mu(T) for any s != 0, and
   * mu -> +inf as det(T) -> 0+, so it forbids element collapse/inversion when
   * driven by Newton from a valid (det > 0) configuration. Zero on rotations
   * and uniform scalings, positive otherwise. value() keeps the simple
   * large-penalty guard for diagnostics; gradient()/hessianAction() implement
   * the exact analytic derivatives on the admissible branch det > 0 (the
   * region the Newton path stays in), so this metric is now usable in the
   * production residual/tangent like SquaredDistanceMetric.
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

      Matrix2 gradient(const Matrix2& T) const
      {
        const Real d = T.determinant();
        const Real F = T.squaredNorm();
        const Matrix2 C = detGradient(T);
        // dmu/dT = T/d - (F / (2 d^2)) * C.
        return (Real(1) / d) * T - (F / (Real(2) * d * d)) * C;
      }

      Matrix2 hessianAction(const Matrix2& T, const Matrix2& H) const
      {
        const Real d = T.determinant();
        const Real F = T.squaredNorm();
        const Matrix2 C = detGradient(T);
        const Matrix2 CofH = detGradient(H);
        const Real ChH = frob(C, H);  // d(det) in direction H
        const Real ThH = frob(T, H);  // (1/2) d(F) in direction H
        const Real invd = Real(1) / d;
        const Real invd2 = invd * invd;
        const Real invd3 = invd2 * invd;
        return invd * H
             - (ChH * invd2) * T
             - (ThH * invd2 - F * ChH * invd3) * C
             - (Real(0.5) * F * invd2) * CofH;
      }

    private:
      // d(det(M))/dM for 2x2 (cofactor matrix); also the linear cofactor map.
      static Matrix2 detGradient(const Matrix2& M)
      {
        Matrix2 C;
        C(0, 0) =  M(1, 1);
        C(0, 1) = -M(1, 0);
        C(1, 0) = -M(0, 1);
        C(1, 1) =  M(0, 0);
        return C;
      }

      static Real frob(const Matrix2& a, const Matrix2& b)
      {
        return a(0, 0) * b(0, 0) + a(0, 1) * b(0, 1)
             + a(1, 0) * b(1, 0) + a(1, 1) * b(1, 1);
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
