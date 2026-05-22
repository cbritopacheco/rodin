/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_METRICS_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_METRICS_H

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/PointCloud.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/H1/H1Element.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  /// Frobenius inner product, any (matching) dimension.
  inline Real matrixInner2D(
      const Math::SpatialMatrix<Real>& a,
      const Math::SpatialMatrix<Real>& b)
  {
    Real s = 0;
    for (std::uint8_t i = 0; i < a.rows(); ++i)
      for (std::uint8_t j = 0; j < a.cols(); ++j)
        s += a(i, j) * b(i, j);
    return s;
  }

  /**
   * @brief Cofactor (derivative of det) of a square matrix, n in {1,2,3}.
   *
   * Returns @f$C@f$ with @f$d\,\det(M)[H] = C:H@f$, i.e. the cofactor matrix
   * (@f$C_{ij} = \partial\det/\partial M_{ij}@f$). Sized from M so the strict
   * TMOP path is dimension-generic; the 2x2 branch is byte-identical to the
   * previous fixed implementation.
   */
  inline Math::SpatialMatrix<Real> cofactor2D(const Math::SpatialMatrix<Real>& M)
  {
    const std::uint8_t n = M.rows();
    Math::SpatialMatrix<Real> C(n, n);
    if (n == 1)
    {
      C(0, 0) = Real(1);
      return C;
    }
    if (n == 2)
    {
      C(0, 0) =  M(1, 1);
      C(0, 1) = -M(1, 0);
      C(1, 0) = -M(0, 1);
      C(1, 1) =  M(0, 0);
      return C;
    }
    for (std::uint8_t i = 0; i < 3; ++i)
      for (std::uint8_t j = 0; j < 3; ++j)
      {
        const std::uint8_t i1 = (i + 1) % 3, i2 = (i + 2) % 3;
        const std::uint8_t j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        C(i, j) = M(i1, j1) * M(i2, j2) - M(i1, j2) * M(i2, j1);
      }
    return C;
  }

  struct Target
  {
    Math::SpatialMatrix<Real> W = Math::SpatialMatrix<Real>::Identity(2, 2);
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
      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto d = static_cast<std::uint8_t>(cell.getDimension());
        return Math::SpatialMatrix<Real>::Identity(d, d);
      }
  };

  inline Math::SpatialMatrix<Real> linearCellJacobian2D(const Geometry::Polytope& cell)
  {
    const auto& mesh = cell.getMesh();
    const auto& vertices = cell.getVertices();
    const auto x0 = mesh.getVertexCoordinates(vertices(0));
    const auto x1 = mesh.getVertexCoordinates(vertices(1));
    const auto x2 = mesh.getVertexCoordinates(vertices(2));

    Math::SpatialMatrix<Real> W(2, 2);
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
        m_targets.resize(mesh.getCellCount(), Math::SpatialMatrix<Real>::Identity(static_cast<std::uint8_t>(mesh.getDimension()), static_cast<std::uint8_t>(mesh.getDimension())));
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

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return linearCellJacobian2D(cell);
      }

    private:
      std::vector<Math::SpatialMatrix<Real>> m_targets;
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
        m_targets.resize(mesh.getCellCount(), Math::SpatialMatrix<Real>::Identity(static_cast<std::uint8_t>(mesh.getDimension()), static_cast<std::uint8_t>(mesh.getDimension())));
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

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return equilateralSameArea(linearCellJacobian2D(cell));
      }

      static Math::SpatialMatrix<Real> equilateralSameArea(const Math::SpatialMatrix<Real>& A0)
      {
        const Real area = Real(0.5) * std::abs(A0.determinant());
        // Degenerate initial element: keep its (tiny) Jacobian so W stays in
        // physical scale instead of collapsing to a singular target.
        if (area <= Real(1e-30))
          return A0;
        const Real l = std::sqrt(Real(4) * area / std::sqrt(Real(3)));
        Math::SpatialMatrix<Real> W(2, 2);
        W(0, 0) = l;
        W(0, 1) = Real(0.5) * l;
        W(1, 0) = Real(0);
        W(1, 1) = l * std::sqrt(Real(3)) / Real(2);
        return W;
      }

    private:
      std::vector<Math::SpatialMatrix<Real>> m_targets;
  };

  /**
   * @brief Per-cell equilateral target with the closest initial orientation.
   *
   * This target keeps the same area prescription as
   * EquilateralSameAreaTargetJacobian, but rotates the equilateral target to
   * be the closest Frobenius-norm fit to the initial element Jacobian. It is
   * better suited to free-boundary moving-mesh solves because neighboring
   * cells are not all driven toward one global reference orientation.
   */
  class OrientedEquilateralSameAreaTargetJacobian
  {
    public:
      OrientedEquilateralSameAreaTargetJacobian() = default;

      template <class Mesh>
      explicit OrientedEquilateralSameAreaTargetJacobian(const Mesh& mesh)
      {
        const auto& conn = mesh.getConnectivity();
        m_targets.resize(mesh.getCellCount(), Math::SpatialMatrix<Real>::Identity(static_cast<std::uint8_t>(mesh.getDimension()), static_cast<std::uint8_t>(mesh.getDimension())));
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          if (conn.getGeometry(2, cellIndex)
              != Geometry::Polytope::Type::Triangle)
            continue;
          auto cellIterator = mesh.getPolytope(2, cellIndex);
          m_targets[static_cast<size_t>(cellIndex)] =
            orientedEquilateralSameArea(linearCellJacobian2D(*cellIterator));
        }
      }

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return orientedEquilateralSameArea(linearCellJacobian2D(cell));
      }

      static Math::SpatialMatrix<Real> orientedEquilateralSameArea(const Math::SpatialMatrix<Real>& A0)
      {
        const Math::SpatialMatrix<Real> E =
          EquilateralSameAreaTargetJacobian::equilateralSameArea(A0);
        if (std::abs(E.determinant()) <= Real(1e-30))
          return E;

        const Math::SpatialMatrix<Real> M = A0 * E.transpose();
        const Real c0 = M(0, 0) + M(1, 1);
        const Real s0 = M(1, 0) - M(0, 1);
        const Real n = std::hypot(c0, s0);
        if (n <= Real(1e-30))
          return E;

        Math::SpatialMatrix<Real> R(2, 2);
        const Real c = c0 / n;
        const Real s = s0 / n;
        R(0, 0) = c;
        R(0, 1) = -s;
        R(1, 0) = s;
        R(1, 1) = c;
        return R * E;
      }

    private:
      std::vector<Math::SpatialMatrix<Real>> m_targets;
  };

  inline Math::SpatialMatrix<Real> parametricCellJacobian2D(
      const Geometry::Polytope& cell);

  /**
   * @brief Closest orientation-preserving regular target with matching measure.
   *
   * This is the quality-improving production target-map abstraction. It does
   * not preserve the current element shape. Instead, each cell receives an
   * ideal same-measure target:
   *
   * - triangles: equilateral same area, rotated to the closest initial
   *   orientation;
   * - quadrilaterals: square same area, rotated to the closest initial
   *   orientation.
   *
   * The current strict TMOP assembly path supports 2D triangles, so the
   * triangle branch is the one exercised today. The quadrilateral branch keeps
   * the target abstraction polytope-oriented for the future without changing
   * the present assembly contract.
   */
  class IdealElementTargetJacobian
  {
    public:
      IdealElementTargetJacobian() = default;

      template <class Mesh>
      explicit IdealElementTargetJacobian(const Mesh& mesh)
      {
        const auto& conn = mesh.getConnectivity();
        const auto d = static_cast<std::uint8_t>(mesh.getDimension());
        m_targets.resize(
            mesh.getCellCount(),
            Math::SpatialMatrix<Real>::Identity(d, d));
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          const auto type = conn.getGeometry(2, cellIndex);
          if (type != Geometry::Polytope::Type::Triangle
              && type != Geometry::Polytope::Type::Quadrilateral)
            continue;
          auto cellIterator = mesh.getPolytope(2, cellIndex);
          m_targets[static_cast<size_t>(cellIndex)] =
            idealSameMeasure(*cellIterator);
        }
      }

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return idealSameMeasure(cell);
      }

      static Math::SpatialMatrix<Real> idealSameMeasure(
          const Geometry::Polytope& cell)
      {
        const auto type =
          cell.getMesh().getGeometry(cell.getDimension(), cell.getIndex());
        const Math::SpatialMatrix<Real> A0 =
          type == Geometry::Polytope::Type::Triangle
            ? linearCellJacobian2D(cell)
            : parametricCellJacobian2D(cell);
        if (type == Geometry::Polytope::Type::Quadrilateral)
          return orientedSquareSameArea(A0);
        return OrientedEquilateralSameAreaTargetJacobian::
          orientedEquilateralSameArea(A0);
      }

      static Math::SpatialMatrix<Real> orientedSquareSameArea(
          const Math::SpatialMatrix<Real>& A0)
      {
        const Real area = std::abs(A0.determinant());
        if (area <= Real(1e-30))
          return A0;

        Math::SpatialMatrix<Real> S(2, 2);
        const Real l = std::sqrt(area);
        S(0, 0) = l;
        S(0, 1) = Real(0);
        S(1, 0) = Real(0);
        S(1, 1) = l;

        const Math::SpatialMatrix<Real> M = A0 * S.transpose();
        const Real c0 = M(0, 0) + M(1, 1);
        const Real s0 = M(1, 0) - M(0, 1);
        const Real n = std::hypot(c0, s0);
        if (n <= Real(1e-30))
          return S;

        Math::SpatialMatrix<Real> R(2, 2);
        const Real c = c0 / n;
        const Real s = s0 / n;
        R(0, 0) = c;
        R(0, 1) = -s;
        R(1, 0) = s;
        R(1, 1) = c;
        return R * S;
      }

    private:
      std::vector<Math::SpatialMatrix<Real>> m_targets;
  };

  /**
   * @brief Generic, polytope- and FES-independent reference Jacobian helper.
   *
   * Uses the cell's parametric map evaluated at the reference-element
   * centroid, so it works for any 2D polytope (triangle, quadrilateral, ...)
   * and any element order without baking in P1/triangle vertex indexing.
   */
  inline Math::SpatialMatrix<Real> parametricCellJacobian2D(const Geometry::Polytope& cell)
  {
    const Geometry::Polytope::Traits traits(cell.getGeometry());
    const Math::SpatialPoint rc = traits.getCentroid();
    const Geometry::Point point(cell, rc);
    const auto& J = point.getJacobian();
    Math::SpatialMatrix<Real> W(2, 2);
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
        : m_W(Math::SpatialMatrix<Real>::Identity(2, 2))
      {}

      explicit ConstantTargetJacobian(const Math::SpatialMatrix<Real>& W)
        : m_W(W)
      {}

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope&,
          const Math::SpatialPoint&) const
      {
        return m_W;
      }

      const Math::SpatialMatrix<Real>& getMatrix() const
      {
        return m_W;
      }

      static ConstantTargetJacobian identity()
      {
        return ConstantTargetJacobian(Math::SpatialMatrix<Real>::Identity(2, 2));
      }

      static ConstantTargetJacobian uniformScale(Real s)
      {
        Math::SpatialMatrix<Real> W = Math::SpatialMatrix<Real>::Identity(2, 2);
        W(0, 0) = s;
        W(1, 1) = s;
        return ConstantTargetJacobian(W);
      }

      static ConstantTargetJacobian diagonal(Real sx, Real sy)
      {
        Math::SpatialMatrix<Real> W = Math::SpatialMatrix<Real>(2, 2);
        W(0, 0) = sx;
        W(1, 1) = sy;
        return ConstantTargetJacobian(W);
      }

      static ConstantTargetJacobian rotation(Real theta)
      {
        const Real c = std::cos(theta);
        const Real s = std::sin(theta);
        Math::SpatialMatrix<Real> W(2, 2);
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
        Math::SpatialMatrix<Real> W(2, 2);
        W(0, 0) = c * sx;
        W(0, 1) = -s * sy;
        W(1, 0) = s * sx;
        W(1, 1) = c * sy;
        return ConstantTargetJacobian(W);
      }

    private:
      Math::SpatialMatrix<Real> m_W;
  };

  /**
   * @brief Spatially varying target driven by a user callable
   * @f$W = f(x)@f$ of the physical point.
   *
   * Polytope- and FES-independent: it only maps the reference sample point to
   * physical coordinates through the cell's parametric map and evaluates the
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

      Math::SpatialMatrix<Real> evaluate(
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
   * Captures, per cell, the initial parametric Jacobian evaluated at the
   * reference centroid via Geometry::Point, so it is correct for triangles,
   * quadrilaterals and higher-order elements alike. On affine triangles it
   * reproduces InitialElementTargetJacobian exactly (and thus gives zero strict
   * energy at @f$u = 0@f$ on the captured mesh), but it does not assume P1
   * triangle vertex indexing.
   */
  class ParametricTargetJacobian
  {
    public:
      ParametricTargetJacobian() = default;

      template <class Mesh>
      explicit ParametricTargetJacobian(const Mesh& mesh)
      {
        const auto& conn = mesh.getConnectivity();
        m_targets.resize(mesh.getCellCount(), Math::SpatialMatrix<Real>::Identity(static_cast<std::uint8_t>(mesh.getDimension()), static_cast<std::uint8_t>(mesh.getDimension())));
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          if (conn.getGeometry(2, cellIndex) == Geometry::Polytope::Type::Point)
            continue;
          auto cellIterator = mesh.getPolytope(2, cellIndex);
          m_targets[static_cast<size_t>(cellIndex)] =
            parametricCellJacobian2D(*cellIterator);
        }
      }

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint&) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (index < m_targets.size())
          return m_targets[index];
        return parametricCellJacobian2D(cell);
      }

    private:
      std::vector<Math::SpatialMatrix<Real>> m_targets;
  };

  /**
   * @brief Quadrature-point target that preserves the current mesh geometry.
   *
   * This is a diagnostic/elastic target, not the production quality-improving
   * target. It evaluates @f$W@f$ as the Jacobian of the current mesh coordinate
   * map at the same reference sample point used by TMOP. Therefore, at zero
   * displacement @f$A = W@f$ and @f$T = A W^{-1} = I@f$ at every sample point.
   * The quality term then regularizes later interface-fit/boundary movement
   * around the current shape instead of driving bad cells toward ideal shapes.
   *
   * Unlike triangle-only targets, this does not hard-code P1 triangle vertex
   * indexing. It uses Geometry::Point and is therefore the natural target-map
   * abstraction for any polytope supported by Rodin's geometry mapping.
   */
  class QualityPreservingTargetJacobian
  {
    public:
      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint& rc) const
      {
        const Geometry::Point point(cell, rc);
        const auto& J = point.getJacobian();
        const auto rows = static_cast<std::uint8_t>(J.rows());
        const auto cols = static_cast<std::uint8_t>(J.cols());
        Math::SpatialMatrix<Real> W(rows, cols);
        for (std::uint8_t i = 0; i < rows; ++i)
          for (std::uint8_t j = 0; j < cols; ++j)
            W(i, j) = J(i, j);
        return W;
      }
  };

  /**
   * @brief Fixed P2 target obtained by projecting interface midside target
   *        nodes onto an analytic feature.
   *
   * This is the curved counterpart to an initial-element target. It does not
   * modify the mesh geometry. Instead, it captures a desired P2 coordinate map
   * at construction time: every cell starts from its affine P2 lift, and only
   * midside nodes on attribute-tagged interface edges are projected by the
   * caller-supplied manifold projector. The quality term then compares
   * @f$A(u)@f$ against a target @f$W@f$ that already contains the intended
   * curvature, so the quality metric no longer asks TMOP to erase the
   * interface curvature that the fit term is trying to recover.
   *
   * Current implementation supports 2D triangular P2 target sampling. Other
   * cells fall back to ParametricTargetJacobian.
   */
  class ProjectedInterfaceTargetJacobian
  {
    public:
      ProjectedInterfaceTargetJacobian() = default;

      template <class Mesh, class Projector>
      ProjectedInterfaceTargetJacobian(
          const Mesh& mesh,
          Geometry::Attribute interfaceAttribute,
          Projector project)
        : m_fallback(mesh),
          m_fe(Geometry::Polytope::Type::Triangle)
      {
        const auto& conn = mesh.getConnectivity();
        m_targets.resize(mesh.getCellCount());
        for (Index ci = 0; ci < static_cast<Index>(mesh.getCellCount()); ++ci)
        {
          if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
            continue;

          const auto affine = affinePointCloud(mesh, ci);
          auto desired = affine;
          for (Index edgeIndex : conn.getIncidence({ 2, 1 }, ci))
          {
            const auto attr = mesh.getAttribute(1, edgeIndex);
            if (!attr || *attr != interfaceAttribute)
              continue;
            const auto local =
              localEdgeVertices(
                  conn.getPolytope(2, ci),
                  conn.getPolytope(1, edgeIndex));
            const size_t node = scalarNodeNear(
                triangleEdgePoint(local[0], local[1], Real(0.5)));
            Math::SpatialPoint x{ affine(0, node), affine(1, node) };
            const auto xp = project(x);
            desired(0, node) = xp[0];
            desired(1, node) = xp[1];
          }

          Real alpha = Real(1);
          for (int attempt = 0; attempt < 20; ++attempt)
          {
            Geometry::PointCloud blended(2, 6);
            for (size_t j = 0; j < 6; ++j)
              for (size_t c = 0; c < 2; ++c)
                blended(c, j) =
                  affine(c, j) + alpha * (desired(c, j) - affine(c, j));
            if (isValid(blended))
            {
              m_targets[static_cast<size_t>(ci)] = std::move(blended);
              break;
            }
            alpha *= Real(0.5);
          }
        }
      }

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint& rc) const
      {
        const auto index = static_cast<size_t>(cell.getIndex());
        if (cell.getGeometry() == Geometry::Polytope::Type::Triangle
            && index < m_targets.size()
            && m_targets[index].getCount() > 0)
        {
          Geometry::ParametricTransformation<
            Variational::RealH1Element<2>> trans(
                m_targets[index], m_fe);
          Math::SpatialMatrix<Real> W;
          trans.jacobian(W, rc);
          if (W.rows() == 2 && W.cols() == 2)
          {
            const Real det = W.determinant();
            if (std::isfinite(det) && det > Real(1e-30))
              return W;
          }
        }
        return m_fallback.evaluate(cell, rc);
      }

    private:
      template <class Mesh>
      Geometry::PointCloud affinePointCloud(const Mesh& mesh, Index ci) const
      {
        const auto& conn = mesh.getConnectivity();
        const auto& ref =
          Variational::RealH1Element<2>::getNodes(
              Geometry::Polytope::Type::Triangle);
        const auto& cell = conn.getPolytope(2, ci);
        const auto a = mesh.getVertexCoordinates(cell(0));
        const auto b = mesh.getVertexCoordinates(cell(1));
        const auto c = mesh.getVertexCoordinates(cell(2));
        Geometry::PointCloud pc(2, ref.size());
        for (size_t j = 0; j < ref.size(); ++j)
        {
          const Real r = ref[j][0];
          const Real s = ref[j][1];
          pc(0, j) = (Real(1) - r - s) * a[0] + r * b[0] + s * c[0];
          pc(1, j) = (Real(1) - r - s) * a[1] + r * b[1] + s * c[1];
        }
        return pc;
      }

      template <class CellKey, class EdgeKey>
      static std::array<std::uint8_t, 2> localEdgeVertices(
          const CellKey& cell,
          const EdgeKey& edge)
      {
        for (std::uint8_t a = 0; a < 3; ++a)
          for (std::uint8_t b = a + 1; b < 3; ++b)
            if ((cell(a) == edge(0) && cell(b) == edge(1))
                || (cell(a) == edge(1) && cell(b) == edge(0)))
              return {{ a, b }};
        return {{ 0, 1 }};
      }

      static Math::SpatialPoint triangleEdgePoint(
          std::uint8_t a,
          std::uint8_t b,
          Real s)
      {
        const auto ca = triangleCorner(a);
        const auto cb = triangleCorner(b);
        return Math::SpatialPoint{
          (Real(1) - s) * ca[0] + s * cb[0],
          (Real(1) - s) * ca[1] + s * cb[1] };
      }

      static Math::SpatialPoint triangleCorner(std::uint8_t k)
      {
        if (k == 1) return Math::SpatialPoint{ Real(1), Real(0) };
        if (k == 2) return Math::SpatialPoint{ Real(0), Real(1) };
        return Math::SpatialPoint{ Real(0), Real(0) };
      }

      static size_t scalarNodeNear(const Math::SpatialPoint& rc)
      {
        const auto& ref =
          Variational::RealH1Element<2>::getNodes(
              Geometry::Polytope::Type::Triangle);
        size_t best = 0;
        Real bestD = std::numeric_limits<Real>::infinity();
        for (size_t j = 0; j < ref.size(); ++j)
        {
          const Real dr = ref[j][0] - rc[0];
          const Real ds = ref[j][1] - rc[1];
          const Real d = dr * dr + ds * ds;
          if (d < bestD)
          {
            bestD = d;
            best = j;
          }
        }
        return best;
      }

      bool isValid(const Geometry::PointCloud& pc) const
      {
        Geometry::ParametricTransformation<
          Variational::RealH1Element<2>> trans(pc, m_fe);
        for (const Math::SpatialPoint& rc : {
              Math::SpatialPoint{ Real(1) / Real(3), Real(1) / Real(3) },
              Math::SpatialPoint{ Real(0.5), Real(0.25) },
              Math::SpatialPoint{ Real(0.25), Real(0.5) },
              Math::SpatialPoint{ Real(0.25), Real(0.25) }})
        {
          Math::SpatialMatrix<Real> J;
          trans.jacobian(J, rc);
          if (J.rows() != 2 || J.cols() != 2)
            return false;
          const Real det = J.determinant();
          if (!std::isfinite(det) || !(det > Real(1e-14)))
            return false;
        }
        return true;
      }

      ParametricTargetJacobian m_fallback;
      Variational::RealH1Element<2> m_fe{
        Geometry::Polytope::Type::Triangle };
      std::vector<Geometry::PointCloud> m_targets;
  };

  /**
   * @brief Curved-geometry target with a quality-improving bias.
   *
   * Curved elements should not be optimized against a purely affine target at
   * every sample point: that asks TMOP to erase curvature while the fit term is
   * trying to place the curved edge on the analytic interface. This target uses
   * the current parametric Jacobian at each quadrature point as the natural
   * curved baseline, then blends it toward an ideal same-measure target for the
   * cell:
   *
   *   W(q) = (1 - theta) J_X(q) + theta W_ideal(K).
   *
   * theta = 0 is purely quality-preserving; theta = 1 is the fully ideal
   * target. Small positive values keep the P2 map natural while still giving a
   * strict quality-improving direction. The blend is locally guarded so the
   * returned target remains orientation preserving whenever the current curved
   * geometry is orientation preserving.
   */
  class CurvedQualityTargetJacobian
  {
    public:
      CurvedQualityTargetJacobian() = default;

      template <class Mesh>
      explicit CurvedQualityTargetJacobian(const Mesh& mesh, Real idealWeight = 0.2)
        : m_ideal(mesh),
          m_idealWeight(clampWeight(idealWeight))
      {}

      CurvedQualityTargetJacobian& setIdealWeight(Real idealWeight)
      {
        m_idealWeight = clampWeight(idealWeight);
        return *this;
      }

      Real getIdealWeight() const
      {
        return m_idealWeight;
      }

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint& rc) const
      {
        const Math::SpatialMatrix<Real> preserve =
          m_preserve.evaluate(cell, rc);
        const Math::SpatialMatrix<Real> ideal =
          m_ideal.evaluate(cell, rc);

        Real theta = m_idealWeight;
        for (int attempt = 0; attempt < 12; ++attempt)
        {
          const Math::SpatialMatrix<Real> W =
            (Real(1) - theta) * preserve + theta * ideal;
          const Real det = W.determinant();
          if (std::isfinite(det) && det > Real(1e-30))
            return W;
          theta *= Real(0.5);
        }
        return preserve;
      }

    private:
      static Real clampWeight(Real idealWeight)
      {
        if (!std::isfinite(idealWeight))
          return Real(0);
        if (idealWeight < Real(0))
          return Real(0);
        if (idealWeight > Real(1))
          return Real(1);
        return idealWeight;
      }

      QualityPreservingTargetJacobian m_preserve;
      IdealElementTargetJacobian m_ideal;
      Real m_idealWeight = 0.2;
  };

  /**
   * @brief Fit-compatible curved target with a controlled quality bias.
   *
   * ProjectedInterfaceTargetJacobian is excellent for not fighting the
   * level-set fit, but it largely preserves the linear cut shape. This target
   * keeps that projected P2 map as the baseline and blends its sampled
   * Jacobian toward the ideal same-measure element:
   *
   *   W(q) = (1 - theta) W_projected(q) + theta W_ideal(K).
   *
   * The important difference from CurvedQualityTargetJacobian is that the
   * baseline already contains the projected interface midside nodes, so the
   * quality term remains compatible with the analytic fit penalty. Small
   * theta values provide a quality-improving direction without asking TMOP to
   * erase the fitted curvature. If a blend would make the target invalid, the
   * weight is halved locally and the projected target is used as the fallback.
   */
  class ProjectedQualityTargetJacobian
  {
    public:
      ProjectedQualityTargetJacobian() = default;

      template <class Mesh, class Projector>
      ProjectedQualityTargetJacobian(
          const Mesh& mesh,
          Geometry::Attribute interfaceAttribute,
          Projector project,
          Real idealWeight = Real(0.05))
        : m_projected(mesh, interfaceAttribute, std::move(project)),
          m_ideal(mesh),
          m_idealWeight(clampWeight(idealWeight))
      {}

      ProjectedQualityTargetJacobian& setIdealWeight(Real idealWeight)
      {
        m_idealWeight = clampWeight(idealWeight);
        return *this;
      }

      Real getIdealWeight() const
      {
        return m_idealWeight;
      }

      Math::SpatialMatrix<Real> evaluate(
          const Geometry::Polytope& cell,
          const Math::SpatialPoint& rc) const
      {
        const Math::SpatialMatrix<Real> projected =
          m_projected.evaluate(cell, rc);
        const Math::SpatialMatrix<Real> ideal =
          m_ideal.evaluate(cell, rc);

        Real theta = m_idealWeight;
        for (int attempt = 0; attempt < 12; ++attempt)
        {
          const Math::SpatialMatrix<Real> W =
            (Real(1) - theta) * projected + theta * ideal;
          const Real det = W.determinant();
          if (std::isfinite(det) && det > Real(1e-30))
            return W;
          theta *= Real(0.5);
        }
        return projected;
      }

    private:
      static Real clampWeight(Real idealWeight)
      {
        if (!std::isfinite(idealWeight))
          return Real(0);
        if (idealWeight < Real(0))
          return Real(0);
        if (idealWeight > Real(1))
          return Real(1);
        return idealWeight;
      }

      ProjectedInterfaceTargetJacobian m_projected;
      IdealElementTargetJacobian m_ideal;
      Real m_idealWeight = Real(0.05);
  };

  class MetricBase
  {
    public:
      virtual ~MetricBase() = default;

      virtual Real value(const Math::SpatialMatrix<Real>& A) const = 0;

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
      Real invalidPenalty(const Math::SpatialMatrix<Real>& A) const
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
      Real value(const Math::SpatialMatrix<Real>& T) const override
      {
        return Real(0.5)
          * (T - Math::SpatialMatrix<Real>::Identity(T.rows(), T.cols())).squaredNorm();
      }

      Math::SpatialMatrix<Real> gradient(const Math::SpatialMatrix<Real>& T) const
      {
        return T - Math::SpatialMatrix<Real>::Identity(T.rows(), T.cols());
      }

      Math::SpatialMatrix<Real> hessianAction(const Math::SpatialMatrix<Real>&, const Math::SpatialMatrix<Real>& dT) const
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
      Real value(const Math::SpatialMatrix<Real>& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        return A.squaredNorm() / (Real(2) * A.determinant()) - Real(1);
      }

      Math::SpatialMatrix<Real> gradient(const Math::SpatialMatrix<Real>& T) const
      {
        const Real d = T.determinant();
        const Real F = T.squaredNorm();
        const Math::SpatialMatrix<Real> C = cofactor2D(T);
        // dmu/dT = T/d - (F / (2 d^2)) * C.
        return (Real(1) / d) * T - (F / (Real(2) * d * d)) * C;
      }

      Math::SpatialMatrix<Real> hessianAction(const Math::SpatialMatrix<Real>& T, const Math::SpatialMatrix<Real>& H) const
      {
        const Real d = T.determinant();
        const Real F = T.squaredNorm();
        const Math::SpatialMatrix<Real> C = cofactor2D(T);
        const Math::SpatialMatrix<Real> CofH = cofactor2D(H);
        const Real ChH = matrixInner2D(C, H);  // d(det) in direction H
        const Real ThH = matrixInner2D(T, H);  // (1/2) d(F) in direction H
        const Real invd = Real(1) / d;
        const Real invd2 = invd * invd;
        const Real invd3 = invd2 * invd;
        return invd * H
             - (ChH * invd2) * T
             - (ThH * invd2 - F * ChH * invd3) * C
             - (Real(0.5) * F * invd2) * CofH;
      }
  };

  /**
   * Area distortion metric:
   *   W(A) = 1/2 (det(A) - 1)^2
   */
  class AreaDistortionMetric final : public MetricBase
  {
    public:
      Real value(const Math::SpatialMatrix<Real>& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        const Real d = A.determinant() - Real(1);
        return Real(0.5) * d * d;
      }

      Math::SpatialMatrix<Real> gradient(const Math::SpatialMatrix<Real>& T) const
      {
        const Real d = T.determinant() - Real(1);
        return d * cofactor2D(T);
      }

      Math::SpatialMatrix<Real> hessianAction(const Math::SpatialMatrix<Real>& T, const Math::SpatialMatrix<Real>& H) const
      {
        const Math::SpatialMatrix<Real> C = cofactor2D(T);
        const Real ddet = matrixInner2D(C, H);
        return ddet * C + (T.determinant() - Real(1)) * cofactor2D(H);
      }
  };

  class ShapeSizeMetric final : public MetricBase
  {
    public:
      Real value(const Math::SpatialMatrix<Real>& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        return m_shape.value(A) + m_area.value(A);
      }

      Math::SpatialMatrix<Real> gradient(const Math::SpatialMatrix<Real>& T) const
      {
        return m_shape.gradient(T) + m_area.gradient(T);
      }

      Math::SpatialMatrix<Real> hessianAction(const Math::SpatialMatrix<Real>& T, const Math::SpatialMatrix<Real>& H) const
      {
        return m_shape.hessianAction(T, H)
             + m_area.hessianAction(T, H);
      }

    private:
      ShapeDistortionMetric m_shape;
      AreaDistortionMetric m_area;
  };

  /**
   * 2D size metric (MFEM TMOP metric 77):
   *   mu(T) = 1/2 (tau - 1/tau)^2,   tau = det(T).
   *
   * Scale-AWARE size barrier: zero iff det(T)=1 (element matches the target
   * size), and mu -> +inf as tau -> 0+, so it forbids size collapse from a
   * valid configuration. This is the size companion used by the
   * shape+size metric mu_80 of Knupp-Kolev-Mittal-Tomov (2021); preferred
   * over (tau-1)^2 because it is symmetric under tau <-> 1/tau and barriers
   * at collapse.
   */
  class SizeMetric77 final : public MetricBase
  {
    public:
      Real value(const Math::SpatialMatrix<Real>& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        const Real tau = A.determinant();
        const Real d = tau - Real(1) / tau;
        return Real(0.5) * d * d;
      }

      Math::SpatialMatrix<Real> gradient(const Math::SpatialMatrix<Real>& T) const
      {
        const Real tau = T.determinant();
        const Real gp = gPrime(tau);
        return gp * cofactor2D(T);
      }

      Math::SpatialMatrix<Real> hessianAction(
          const Math::SpatialMatrix<Real>& T,
          const Math::SpatialMatrix<Real>& H) const
      {
        const Real tau = T.determinant();
        const Math::SpatialMatrix<Real> C = cofactor2D(T);
        const Real dtau = matrixInner2D(C, H);    // d(det) in direction H
        return gPrimePrime(tau) * dtau * C + gPrime(tau) * cofactor2D(H);
      }

    private:
      // g(tau) = 1/2 (tau - 1/tau)^2.
      static Real gPrime(Real tau)
      {
        return (tau - Real(1) / tau) * (Real(1) + Real(1) / (tau * tau));
      }
      static Real gPrimePrime(Real tau)
      {
        const Real inv2 = Real(1) / (tau * tau);
        return (Real(1) + inv2) * (Real(1) + inv2)
             + (tau - Real(1) / tau) * (-Real(2) / (tau * tau * tau));
      }
  };

  /**
   * 2D shape+size metric (MFEM TMOP metric 80):
   *   mu_80 = (1 - gamma) mu_2 + gamma mu_77,   gamma in [0,1].
   *
   * mu_2  = ShapeDistortionMetric (scale-invariant shape barrier),
   * mu_77 = SizeMetric77          (size barrier).
   *
   * This is the production metric for moving-interface fitting in
   * Knupp-Kolev-Mittal-Tomov (2021) (their Taylor-Green / Rayleigh-Taylor
   * tests use gamma=0.5). Pure-shape metrics are scale-rank-deficient and
   * leave slivers; the size term restores control over element size.
   */
  class ShapeSizeBlendMetric final : public MetricBase
  {
    public:
      explicit ShapeSizeBlendMetric(Real gamma = Real(0.5))
        : m_gamma(std::max(Real(0), std::min(Real(1), gamma)))
      {}

      Real getGamma() const { return m_gamma; }

      Real value(const Math::SpatialMatrix<Real>& A) const override
      {
        if (const Real penalty = invalidPenalty(A); penalty > Real(0))
          return penalty;
        return (Real(1) - m_gamma) * m_shape.value(A)
             + m_gamma * m_size.value(A);
      }

      Math::SpatialMatrix<Real> gradient(const Math::SpatialMatrix<Real>& T) const
      {
        return (Real(1) - m_gamma) * m_shape.gradient(T)
             + m_gamma * m_size.gradient(T);
      }

      Math::SpatialMatrix<Real> hessianAction(
          const Math::SpatialMatrix<Real>& T,
          const Math::SpatialMatrix<Real>& H) const
      {
        return (Real(1) - m_gamma) * m_shape.hessianAction(T, H)
             + m_gamma * m_size.hessianAction(T, H);
      }

    private:
      Real m_gamma;
      ShapeDistortionMetric m_shape;  // mu_2
      SizeMetric77 m_size;            // mu_77
  };
}

#endif
