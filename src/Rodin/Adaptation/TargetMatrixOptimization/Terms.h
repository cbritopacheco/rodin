/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_TERMS_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_TERMS_H

/**
 * @file
 * @brief Rodin-native TMOP residual and tangent terms.
 *
 * These terms are intentionally thin adapters over Rodin's variational
 * infrastructure. They return normal linear-form and bilinear-form integrators,
 * so assembly, linear solves, Newton iterations, PETSc, MPI, and OpenMP support
 * keep flowing through the existing LinearProblem machinery.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Rodin/Adaptation/TargetMatrixOptimization/Metrics.h"
#include "Rodin/Geometry/LevelSetDiscretizerTriangles.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  inline Real triangleArea2D(
        const Math::SpatialPoint& x0,
        const Math::SpatialPoint& x1,
        const Math::SpatialPoint& x2)
    {
      return Real(0.5) * std::abs(
          (x1[0] - x0[0]) * (x2[1] - x0[1])
        - (x1[1] - x0[1]) * (x2[0] - x0[0]));
    }

  inline Real edgeLength(
        const Math::SpatialPoint& x0,
        const Math::SpatialPoint& x1)
    {
      return (x1 - x0).norm();
    }

  inline Real equilateralEdgeLengthFromArea(Real area)
    {
      if (area <= Real(0))
        return Real(0);
      return std::sqrt(Real(4) * area / std::sqrt(Real(3)));
    }

  template <class GridFunctionType>
  Math::SpatialPoint deformedVertex(
        const GridFunctionType& u,
        const Geometry::Polytope& polytope,
        size_t localVertex,
        size_t sdim)
    {
      const auto& mesh = u.getFiniteElementSpace().getMesh();
      const Index vertex = polytope.getVertices()(localVertex);
      auto x = mesh.getVertexCoordinates(vertex);
      const auto& data = u.getData();
      const Index vertexCount = static_cast<Index>(mesh.getVertexCount());
      for (size_t c = 0; c < sdim; ++c)
        x[c] += data(vertex + static_cast<Index>(c) * vertexCount);
      return x;
    }

  inline Real frobeniusInner(const Math::SpatialMatrix<Real>& a, const Math::SpatialMatrix<Real>& b)
  {
    Real s = 0;
    for (std::uint8_t i = 0; i < a.rows(); ++i)
      for (std::uint8_t j = 0; j < a.cols(); ++j)
        s += a(i, j) * b(i, j);
    return s;
  }

  // Copies an arbitrary (square) matrix-like into a SpatialMatrix sized from
  // the source, so the strict path is dimension-generic. 2D is unchanged.
  template <class MatrixType>
  Math::SpatialMatrix<Real> toSpatialMatrix(const MatrixType& matrix)
  {
    const auto r = static_cast<std::uint8_t>(matrix.rows());
    const auto c = static_cast<std::uint8_t>(matrix.cols());
    Math::SpatialMatrix<Real> out(r, c);
    for (std::uint8_t i = 0; i < r; ++i)
      for (std::uint8_t j = 0; j < c; ++j)
        out(i, j) = matrix(i, j);
    return out;
  }

  inline void requireStrictTMOPCell(
      const Geometry::Polytope& cell,
      size_t vectorDimension,
      size_t localSize)
  {
    const auto geometry =
      cell.getMesh().getGeometry(cell.getDimension(), cell.getIndex());
    if (cell.getDimension() != 2
        || geometry != Geometry::Polytope::Type::Triangle
        || vectorDimension != 2
        || localSize == 0
        || localSize % vectorDimension != 0)
    {
      throw std::runtime_error(
          "TMOP::QualityTerm currently supports 2D triangular vector "
          "finite element spaces with two displacement components.");
    }
  }

  template <class ElementType>
  Math::SpatialMatrix<Real> basisJacobian2D(
      const ElementType& fe,
      size_t local,
      const Math::SpatialPoint& rc)
  {
    return toSpatialMatrix(fe.getBasis(local).getJacobian()(rc));
  }

  template <class GridFunctionType>
  std::vector<Real> localCoordinateCoefficients(
      const GridFunctionType& u,
      const Geometry::Polytope& cell)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& fe = fes.getFiniteElement(cell.getDimension(), cell.getIndex());
    const size_t vdim = fes.getVectorDimension();
    const size_t localSize = fe.getCount();
    requireStrictTMOPCell(cell, vdim, localSize);

    std::vector<Real> coeff(localSize, Real(0));
    const auto& data = u.getData();
    const size_t scalarLocalSize = localSize / vdim;
    for (size_t node = 0; node < scalarLocalSize; ++node)
    {
      const auto& nodeRef = fe.getNode(node * vdim);
      const Geometry::Point point(cell, nodeRef);
      const auto& X = point.getPhysicalCoordinates();
      for (size_t component = 0; component < vdim; ++component)
      {
        const size_t local = node * vdim + component;
        const Index global = fes.getGlobalIndex(
            {cell.getDimension(), cell.getIndex()},
            static_cast<Index>(local));
        coeff[local] = X[component] + data(global);
      }
    }
    return coeff;
  }

  template <class ElementType>
  Math::SpatialMatrix<Real> deformedCoordinateJacobian(
      const ElementType& fe,
      const std::vector<Real>& coefficients,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    Math::SpatialMatrix<Real> A(static_cast<std::uint8_t>(vdim),
              static_cast<std::uint8_t>(vdim));
    for (size_t local = 0; local < coefficients.size(); ++local)
      A += coefficients[local] * basisJacobian2D(fe, local, rc);
    return A;
  }

  template <class GridFunctionType>
  Math::SpatialMatrix<Real> deformedCoordinateJacobian(
      const GridFunctionType& u,
      const Geometry::Polytope& cell,
      const Math::SpatialPoint& rc)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& fe = fes.getFiniteElement(cell.getDimension(), cell.getIndex());
    return deformedCoordinateJacobian(
        fe, localCoordinateCoefficients(u, cell), rc, fes.getVectorDimension());
  }

  template <class ElementType>
  Real basisComponentValue(
      const ElementType& fe,
      size_t local,
      const Math::SpatialPoint& rc,
      size_t component)
  {
    return fe.getBasis(local)(rc)[static_cast<std::uint8_t>(component)];
  }

  template <class ElementType>
  Real basisComponentDerivative1D(
      const ElementType& fe,
      size_t local,
      const Math::SpatialPoint& rc,
      size_t component)
  {
    const auto J = fe.getBasis(local).getJacobian()(rc);
    return J(static_cast<std::uint8_t>(component), 0);
  }

  template <class GridFunctionType>
  std::vector<Real> localEdgeCoordinateCoefficients(
      const GridFunctionType& u,
      const Geometry::Polytope& edge)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& fe = fes.getFiniteElement(edge.getDimension(), edge.getIndex());
    const size_t vdim = fes.getVectorDimension();
    const size_t localSize = fe.getCount();
    if (edge.getDimension() != 1
        || edge.getMesh().getSpaceDimension() != 2
        || vdim != 2
        || localSize == 0
        || localSize % vdim != 0)
    {
      throw std::runtime_error(
          "TMOP::AnalyticLevelSetFitTerm currently supports 2D vector "
          "finite element spaces on segment interface edges.");
    }

    std::vector<Real> coeff(localSize, Real(0));
    const auto& data = u.getData();
    const size_t scalarLocalSize = localSize / vdim;
    for (size_t node = 0; node < scalarLocalSize; ++node)
    {
      const auto& nodeRef = fe.getNode(node * vdim);
      const Geometry::Point point(edge, nodeRef);
      const auto& X = point.getPhysicalCoordinates();
      for (size_t component = 0; component < vdim; ++component)
      {
        const size_t local = node * vdim + component;
        const Index global = fes.getGlobalIndex(
            {edge.getDimension(), edge.getIndex()},
            static_cast<Index>(local));
        coeff[local] = X[component] + data(global);
      }
    }
    return coeff;
  }

  template <class ElementType>
  Math::SpatialPoint deformedEdgePoint(
      const ElementType& fe,
      const std::vector<Real>& coefficients,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    Math::SpatialPoint x(static_cast<std::uint8_t>(vdim));
    x.setZero();
    for (size_t local = 0; local < coefficients.size(); ++local)
    {
      const size_t component = local % vdim;
      x[static_cast<std::uint8_t>(component)] +=
        coefficients[local] * basisComponentValue(fe, local, rc, component);
    }
    return x;
  }

  template <class ElementType>
  Math::SpatialPoint deformedEdgeDerivative(
      const ElementType& fe,
      const std::vector<Real>& coefficients,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    Math::SpatialPoint dx(static_cast<std::uint8_t>(vdim));
    dx.setZero();
    for (size_t local = 0; local < coefficients.size(); ++local)
    {
      const size_t component = local % vdim;
      dx[static_cast<std::uint8_t>(component)] +=
        coefficients[local]
      * basisComponentDerivative1D(fe, local, rc, component);
    }
    return dx;
  }

  template <class GridFunctionType, class Metric, class TargetJacobian>
  class TMOPQualityResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LinearFormIntegratorBase<Real>;

        template <class TestFES>
        TMOPQualityResidualIntegrator(
            const GridFunctionType& u,
            const Variational::TestFunction<TestFES>& v,
            Metric metric,
            TargetJacobian target,
            Real weight,
            size_t quadratureOrder)
          : Parent(v),
            m_u(u),
            m_metric(std::move(metric)),
            m_target(std::move(target)),
            m_weight(std::max(Real(0), weight)),
            m_quadratureOrder(quadratureOrder)
        {}

        TMOPQualityResidualIntegrator(
            const TMOPQualityResidualIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_metric(other.m_metric),
            m_target(other.m_target),
            m_weight(other.m_weight),
            m_quadratureOrder(other.m_quadratureOrder),
            m_polytope(other.m_polytope)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        // Per-cell cache. The deformed Jacobian, target, T and metric
        // gradient are functions of the cell+quadrature only (NOT of the test
        // dof), and basisJacobian2D depends only on (local, q). Computing them
        // once per cell here turns integrate(local) into a cheap contraction
        // instead of recomputing the heavy Geometry::Point-based Jacobian for
        // every local dof (the dominant assembly cost).
        TMOPQualityResidualIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          const auto& cell = polytope;
          const auto& fes = m_u.get().getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          requireStrictTMOPCell(cell, fes.getVectorDimension(), fe.getCount());
          m_localCount = fe.getCount();
          const auto coefficients = localCoordinateCoefficients(m_u.get(), cell);
          const auto geometry =
            cell.getMesh().getGeometry(cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(m_quadratureOrder, geometry);
          const size_t nq = qf.getSize();
          m_wdetW.assign(nq, 0);
          const auto d = static_cast<std::uint8_t>(fes.getVectorDimension());
          m_G.assign(nq, Math::SpatialMatrix<Real>(d, d));
          m_dT.assign(nq, std::vector<Math::SpatialMatrix<Real>>(m_localCount));
          for (size_t q = 0; q < nq; ++q)
          {
            const auto& rc = qf.getPoint(q);
            std::vector<Math::SpatialMatrix<Real>> basis(m_localCount);
            Math::SpatialMatrix<Real> A(d, d);
            for (size_t l = 0; l < m_localCount; ++l)
            {
              basis[l] = basisJacobian2D(fe, l, rc);
              A += coefficients[l] * basis[l];
            }
            const Math::SpatialMatrix<Real> W = m_target.evaluate(cell, rc);
            const Real detW = W.determinant();
            if (std::abs(detW) <= Real(1e-30))
              throw std::runtime_error(
                  "TMOP::QualityTerm target Jacobian is singular.");
            const Math::SpatialMatrix<Real> Winv = W.inverse();
            m_wdetW[q] = qf.getWeight(q) * detW;
            m_G[q] = m_metric.gradient(A * Winv);
            for (size_t l = 0; l < m_localCount; ++l)
              m_dT[q][l] = basis[l] * Winv;
          }
          return *this;
        }

        Real integrate(size_t local) override
        {
          if (local >= m_localCount)
            return 0;
          Real value = 0;
          for (size_t q = 0; q < m_wdetW.size(); ++q)
            value += m_wdetW[q]
              * frobeniusInner(m_G[q], m_dT[q][local]);
          return m_weight * value;
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Cells;
        }

        TMOPQualityResidualIntegrator* copy() const noexcept override
        {
          return new TMOPQualityResidualIntegrator(*this);
        }

      private:
        std::reference_wrapper<const GridFunctionType> m_u;
        Metric m_metric;
        TargetJacobian m_target;
        Real m_weight = 1;
        size_t m_quadratureOrder = 2;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localCount = 0;
        std::vector<Real> m_wdetW;
        std::vector<Math::SpatialMatrix<Real>> m_G;
        std::vector<std::vector<Math::SpatialMatrix<Real>>> m_dT;
    };

  template <class GridFunctionType, class Metric, class TargetJacobian>
  class TMOPQualityTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

        template <class Solution, class TrialFES, class TestFES>
        TMOPQualityTangentIntegrator(
            const GridFunctionType& u,
            const Variational::TrialFunction<Solution, TrialFES>& du,
            const Variational::TestFunction<TestFES>& v,
            Metric metric,
            TargetJacobian target,
            Real weight,
            size_t quadratureOrder)
          : Parent(du, v),
            m_u(u),
            m_metric(std::move(metric)),
            m_target(std::move(target)),
            m_weight(std::max(Real(0), weight)),
            m_quadratureOrder(quadratureOrder)
        {}

        TMOPQualityTangentIntegrator(
            const TMOPQualityTangentIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_metric(other.m_metric),
            m_target(other.m_target),
            m_weight(other.m_weight),
            m_quadratureOrder(other.m_quadratureOrder),
            m_polytope(other.m_polytope)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        // Per-cell cache: A, T, target and the metric Hessian action on every
        // trial dof are functions of (cell, quadrature) only. Precomputing
        // them here makes integrate(trial,test) an O(nq) contraction instead
        // of recomputing the Geometry::Point-based deformed Jacobian and the
        // metric Hessian for every one of the O(localDofs^2) dof pairs.
        TMOPQualityTangentIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          const auto& cell = polytope;
          const auto& fes = m_u.get().getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          requireStrictTMOPCell(cell, fes.getVectorDimension(), fe.getCount());
          m_localCount = fe.getCount();
          const auto coefficients = localCoordinateCoefficients(m_u.get(), cell);
          const auto geometry =
            cell.getMesh().getGeometry(cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(m_quadratureOrder, geometry);
          const size_t nq = qf.getSize();
          m_wdetW.assign(nq, 0);
          m_dT.assign(nq, std::vector<Math::SpatialMatrix<Real>>(m_localCount));
          m_Hd.assign(nq, std::vector<Math::SpatialMatrix<Real>>(m_localCount));
          for (size_t q = 0; q < nq; ++q)
          {
            const auto& rc = qf.getPoint(q);
            std::vector<Math::SpatialMatrix<Real>> basis(m_localCount);
            Math::SpatialMatrix<Real> A(
                static_cast<std::uint8_t>(fes.getVectorDimension()),
                static_cast<std::uint8_t>(fes.getVectorDimension()));
            for (size_t l = 0; l < m_localCount; ++l)
            {
              basis[l] = basisJacobian2D(fe, l, rc);
              A += coefficients[l] * basis[l];
            }
            const Math::SpatialMatrix<Real> W = m_target.evaluate(cell, rc);
            const Real detW = W.determinant();
            if (std::abs(detW) <= Real(1e-30))
              throw std::runtime_error(
                  "TMOP::QualityTerm target Jacobian is singular.");
            const Math::SpatialMatrix<Real> Winv = W.inverse();
            const Math::SpatialMatrix<Real> T = A * Winv;
            m_wdetW[q] = qf.getWeight(q) * detW;
            for (size_t l = 0; l < m_localCount; ++l)
              m_dT[q][l] = basis[l] * Winv;
            for (size_t l = 0; l < m_localCount; ++l)
              m_Hd[q][l] = m_metric.hessianAction(T, m_dT[q][l]);
          }
          return *this;
        }

        Real integrate(size_t trial, size_t test) override
        {
          if (trial >= m_localCount || test >= m_localCount)
            return 0;
          Real value = 0;
          for (size_t q = 0; q < m_wdetW.size(); ++q)
            value += m_wdetW[q]
              * frobeniusInner(m_Hd[q][trial], m_dT[q][test]);
          return m_weight * value;
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Cells;
        }

        TMOPQualityTangentIntegrator* copy() const noexcept override
        {
          return new TMOPQualityTangentIntegrator(*this);
        }

      private:
        std::reference_wrapper<const GridFunctionType> m_u;
        Metric m_metric;
        TargetJacobian m_target;
        Real m_weight = 1;
        size_t m_quadratureOrder = 2;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localCount = 0;
        std::vector<Real> m_wdetW;
        std::vector<std::vector<Math::SpatialMatrix<Real>>> m_dT;
        std::vector<std::vector<Math::SpatialMatrix<Real>>> m_Hd;
    };

  template <class GridFunctionType>
  class EdgeSpringQualityResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LinearFormIntegratorBase<Real>;

        template <class TestFES>
        EdgeSpringQualityResidualIntegrator(
            const GridFunctionType& u,
            const Variational::TestFunction<TestFES>& v,
            Real weight)
          : Parent(v),
            m_u(u),
            m_weight(std::max(Real(0), weight))
        {}

        EdgeSpringQualityResidualIntegrator(
            const EdgeSpringQualityResidualIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_local(other.m_local)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        EdgeSpringQualityResidualIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_local.clear();

          if (polytope.getDimension() != 2
              || polytope.getVertices().size() != 3)
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          m_local = computeLocalResidual(polytope, sdim);
          return *this;
        }

        Real integrate(size_t local) override
        {
          if (local >= m_local.size())
            return 0;
          return m_local[local];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Cells;
        }

        EdgeSpringQualityResidualIntegrator* copy() const noexcept override
        {
          return new EdgeSpringQualityResidualIntegrator(*this);
        }

      private:
        std::vector<Real> computeLocalResidual(
            const Geometry::Polytope& cell,
            size_t sdim) const
        {
          std::vector<Real> local(cell.getVertices().size() * sdim, Real(0));

          std::array<Math::SpatialPoint, 3> X;
          std::array<Math::SpatialPoint, 3> x;
          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          for (size_t i = 0; i < 3; ++i)
          {
            X[i] = mesh.getVertexCoordinates(cell.getVertices()(i));
            x[i] = deformedVertex(m_u.get(), cell, i, sdim);
          }

          const Real area = triangleArea2D(X[0], X[1], X[2]);
          Real target = equilateralEdgeLengthFromArea(area);
          if (target <= Real(0))
          {
            target =
              (edgeLength(X[0], X[1])
             + edgeLength(X[1], X[2])
             + edgeLength(X[2], X[0])) / Real(3);
          }

          static constexpr std::array<std::array<size_t, 2>, 3> Edges = {{
            {{0, 1}}, {{1, 2}}, {{2, 0}}
          }};

          for (const auto& edge : Edges)
          {
            const size_t a = edge[0];
            const size_t b = edge[1];
            const auto e = x[b] - x[a];
            const Real length = e.norm();
            if (length <= Real(1e-14))
              continue;
            const auto direction = (Real(1) / length) * e;
            const Real scale = m_weight * (length - target);
            for (size_t c = 0; c < sdim; ++c)
            {
              local[a * sdim + c] -= scale * direction[c];
              local[b * sdim + c] += scale * direction[c];
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        std::vector<Real> m_local;
    };

  template <class GridFunctionType>
  class EdgeSpringQualityTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

        template <class Solution, class TrialFES, class TestFES>
        EdgeSpringQualityTangentIntegrator(
            const GridFunctionType& u,
            const Variational::TrialFunction<Solution, TrialFES>& du,
            const Variational::TestFunction<TestFES>& v,
            Real weight)
          : Parent(du, v),
            m_u(u),
            m_weight(std::max(Real(0), weight))
        {}

        EdgeSpringQualityTangentIntegrator(
            const EdgeSpringQualityTangentIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_localSize(other.m_localSize),
            m_matrix(other.m_matrix)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        EdgeSpringQualityTangentIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_localSize = 0;
          m_matrix.clear();

          if (polytope.getDimension() != 2
              || polytope.getVertices().size() != 3)
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          m_localSize = polytope.getVertices().size() * sdim;
          m_matrix = computeLocalTangent(polytope, sdim);
          return *this;
        }

        Real integrate(size_t trial, size_t test) override
        {
          if (trial >= m_localSize || test >= m_localSize)
            return 0;
          return m_matrix[test * m_localSize + trial];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Cells;
        }

        EdgeSpringQualityTangentIntegrator* copy() const noexcept override
        {
          return new EdgeSpringQualityTangentIntegrator(*this);
        }

      private:
        std::vector<Real> computeLocalTangent(
            const Geometry::Polytope& cell,
            size_t sdim) const
        {
          const size_t localSize = cell.getVertices().size() * sdim;
          std::vector<Real> local(localSize * localSize, Real(0));

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          std::array<Math::SpatialPoint, 3> X;
          std::array<Math::SpatialPoint, 3> x;
          for (size_t i = 0; i < 3; ++i)
          {
            X[i] = mesh.getVertexCoordinates(cell.getVertices()(i));
            x[i] = deformedVertex(m_u.get(), cell, i, sdim);
          }

          // Same rest length as the residual (equilateral length from the
          // undeformed area, with a fallback), so the tangent below is the
          // exact Jacobian of the residual edge force.
          const Real area = triangleArea2D(X[0], X[1], X[2]);
          Real target = equilateralEdgeLengthFromArea(area);
          if (target <= Real(0))
          {
            target =
              (edgeLength(X[0], X[1])
             + edgeLength(X[1], X[2])
             + edgeLength(X[2], X[0])) / Real(3);
          }

          static constexpr std::array<std::array<size_t, 2>, 3> Edges = {{
            {{0, 1}}, {{1, 2}}, {{2, 0}}
          }};

          auto add = [&](size_t rowVertex, size_t rowComp,
                         size_t colVertex, size_t colComp,
                         Real value)
          {
            const size_t row = rowVertex * sdim + rowComp;
            const size_t col = colVertex * sdim + colComp;
            local[row * localSize + col] += value;
          };

          for (const auto& edge : Edges)
          {
            const size_t a = edge[0];
            const size_t b = edge[1];
            const auto e = x[b] - x[a];
            const Real length = e.norm();
            if (length <= Real(1e-14))
              continue;
            const auto direction = (Real(1) / length) * e;
            // d/dx of f = w (|e| - L0) e/|e| :
            //   K = w [ ee^T/|e|^2  +  (|e|-L0)/|e| (I - ee^T/|e|^2) ].
            // The first term is the axial stiffness; the second is the
            // geometric stiffness that the previous tangent omitted, which
            // is what broke residual/Jacobian consistency away from
            // equilibrium and forced tiny damping.
            const Real geom = (length - target) / length;
            for (size_t r = 0; r < sdim; ++r)
            {
              for (size_t cc = 0; cc < sdim; ++cc)
              {
                const Real dd = direction[r] * direction[cc];
                const Real identity = (r == cc) ? Real(1) : Real(0);
                const Real value =
                  m_weight * (dd + geom * (identity - dd));
                add(a, r, a, cc,  value);
                add(a, r, b, cc, -value);
                add(b, r, a, cc, -value);
                add(b, r, b, cc,  value);
              }
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localSize = 0;
        std::vector<Real> m_matrix;
    };

  template <class GridFunctionType>
  class SourceSegmentFitResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LinearFormIntegratorBase<Real>;

        template <class TestFES>
        SourceSegmentFitResidualIntegrator(
            const GridFunctionType& u,
            const Variational::TestFunction<TestFES>& v,
            const Geometry::InterfaceGraph& graph,
            const Geometry::LevelSetDiscretizerTrianglesReport& report,
            Real weight)
          : Parent(v),
            m_u(u),
            m_graph(graph),
            m_report(report),
            m_weight(std::max(Real(0), weight))
        {}

        SourceSegmentFitResidualIntegrator(
            const SourceSegmentFitResidualIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_graph(other.m_graph),
            m_report(other.m_report),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_local(other.m_local)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        SourceSegmentFitResidualIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_local.clear();

          if (polytope.getDimension() != 1
              || polytope.getVertices().size() != 2)
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          m_local = computeLocalResidual(polytope, sdim);
          return *this;
        }

        Real integrate(size_t local) override
        {
          if (local >= m_local.size())
            return 0;
          return m_local[local];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Faces;
        }

        SourceSegmentFitResidualIntegrator* copy() const noexcept override
        {
          return new SourceSegmentFitResidualIntegrator(*this);
        }

      private:
        std::vector<Real> computeLocalResidual(
            const Geometry::Polytope& edge,
            size_t sdim) const
        {
          std::vector<Real> local(edge.getVertices().size() * sdim, Real(0));
          const auto provenanceIt =
            m_report.get().interfaceEdgeProvenance.find(edge.getIndex());
          if (provenanceIt == m_report.get().interfaceEdgeProvenance.end())
            return local;

          const auto& provenance = provenanceIt->second;
          if (provenance.sourceInterfaceGraphEdge >= m_graph.get().edges.size())
            return local;

          const auto& source =
            m_graph.get().edges[provenance.sourceInterfaceGraphEdge];
          const auto& a = m_graph.get().vertices[source.v0].x;
          const auto& b = m_graph.get().vertices[source.v1].x;
          const auto ab = b - a;
          const Real denom = ab.dot(ab);
          if (denom <= Real(1e-30))
            return local;

          auto x0 = deformedVertex(m_u.get(), edge, 0, sdim);
          auto x1 = deformedVertex(m_u.get(), edge, 1, sdim);
          const Real length = (x1 - x0).norm();
          if (length <= Real(1e-14))
            return local;

          const auto t = (Real(1) / std::sqrt(denom)) * ab;

          static constexpr std::array<Real, 2> Points = {{
            Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
            Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
          }};
          static constexpr Real Weight = Real(0.5);

          for (Real q : Points)
          {
            const std::array<Real, 2> N = {{Real(1) - q, q}};
            const auto x = N[0] * x0 + N[1] * x1;
            auto d = x - a;
            d -= d.dot(t) * t;
            for (size_t node = 0; node < 2; ++node)
              for (size_t c = 0; c < sdim; ++c)
                local[node * sdim + c] +=
                  m_weight * Weight * length * N[node] * d[c];
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        std::reference_wrapper<const Geometry::InterfaceGraph> m_graph;
        std::reference_wrapper<const Geometry::LevelSetDiscretizerTrianglesReport> m_report;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        std::vector<Real> m_local;
    };

  template <class GridFunctionType>
  class DeviationResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LinearFormIntegratorBase<Real>;

        template <class TestFES>
        DeviationResidualIntegrator(
            const GridFunctionType& u,
            const Variational::TestFunction<TestFES>& v,
            Real weight)
          : Parent(v),
            m_u(u),
            m_weight(std::max(Real(0), weight))
        {}

        DeviationResidualIntegrator(
            const DeviationResidualIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_local(other.m_local)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        DeviationResidualIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_local.clear();

          if (polytope.getDimension() != 2
              || polytope.getVertices().size() != 3)
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          m_local = computeLocalResidual(polytope, sdim);
          return *this;
        }

        Real integrate(size_t local) override
        {
          if (local >= m_local.size())
            return 0;
          return m_local[local];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Cells;
        }

        DeviationResidualIntegrator* copy() const noexcept override
        {
          return new DeviationResidualIntegrator(*this);
        }

      private:
        std::vector<Real> computeLocalResidual(
            const Geometry::Polytope& cell,
            size_t sdim) const
        {
          const size_t localSize = cell.getVertices().size() * sdim;
          std::vector<Real> local(localSize, Real(0));

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const auto x0 = mesh.getVertexCoordinates(cell.getVertices()(0));
          const auto x1 = mesh.getVertexCoordinates(cell.getVertices()(1));
          const auto x2 = mesh.getVertexCoordinates(cell.getVertices()(2));
          const Real area = triangleArea2D(x0, x1, x2);
          const auto& data = m_u.get().getData();
          const Index vertexCount = static_cast<Index>(mesh.getVertexCount());

          for (size_t a = 0; a < 3; ++a)
          {
            for (size_t b = 0; b < 3; ++b)
            {
              const Real mass = m_weight * area
                * ((a == b) ? Real(1) / Real(6) : Real(1) / Real(12));
              const Index vertexB = cell.getVertices()(b);
              for (size_t c = 0; c < sdim; ++c)
              {
                local[a * sdim + c] +=
                  mass * data(vertexB + static_cast<Index>(c) * vertexCount);
              }
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        std::vector<Real> m_local;
    };

  class DeviationTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

        template <class Solution, class TrialFES, class TestFES>
        DeviationTangentIntegrator(
            const Variational::TrialFunction<Solution, TrialFES>& du,
            const Variational::TestFunction<TestFES>& v,
            Real weight)
          : Parent(du, v),
            m_weight(std::max(Real(0), weight))
        {}

        DeviationTangentIntegrator(
            const DeviationTangentIntegrator& other)
          : Parent(other),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_localSize(other.m_localSize),
            m_matrix(other.m_matrix)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        DeviationTangentIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_localSize = 0;
          m_matrix.clear();

          if (polytope.getDimension() != 2
              || polytope.getVertices().size() != 3)
            return *this;

          const size_t sdim = 2;
          m_localSize = polytope.getVertices().size() * sdim;
          m_matrix.assign(m_localSize * m_localSize, Real(0));

          const auto& cell = polytope;
          const auto& mesh = cell.getMesh();
          const auto x0 = mesh.getVertexCoordinates(cell.getVertices()(0));
          const auto x1 = mesh.getVertexCoordinates(cell.getVertices()(1));
          const auto x2 = mesh.getVertexCoordinates(cell.getVertices()(2));
          const Real area = triangleArea2D(x0, x1, x2);

          for (size_t trialNode = 0; trialNode < 3; ++trialNode)
          {
            for (size_t testNode = 0; testNode < 3; ++testNode)
            {
              const Real mass = m_weight * area
                * ((trialNode == testNode)
                    ? Real(1) / Real(6)
                    : Real(1) / Real(12));
              for (size_t c = 0; c < sdim; ++c)
              {
                const size_t row = testNode * sdim + c;
                const size_t col = trialNode * sdim + c;
                m_matrix[row * m_localSize + col] = mass;
              }
            }
          }
          return *this;
        }

        Real integrate(size_t trial, size_t test) override
        {
          if (trial >= m_localSize || test >= m_localSize)
            return 0;
          return m_matrix[test * m_localSize + trial];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Cells;
        }

        DeviationTangentIntegrator* copy() const noexcept override
        {
          return new DeviationTangentIntegrator(*this);
        }

      private:
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localSize = 0;
        std::vector<Real> m_matrix;
    };

  template <class GridFunctionType>
  class SourceSegmentFitTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

        template <class Solution, class TrialFES, class TestFES>
        SourceSegmentFitTangentIntegrator(
            const GridFunctionType& u,
            const Variational::TrialFunction<Solution, TrialFES>& du,
            const Variational::TestFunction<TestFES>& v,
            const Geometry::InterfaceGraph& graph,
            const Geometry::LevelSetDiscretizerTrianglesReport& report,
            Real weight)
          : Parent(du, v),
            m_u(u),
            m_graph(graph),
            m_report(report),
            m_weight(std::max(Real(0), weight))
        {}

        SourceSegmentFitTangentIntegrator(
            const SourceSegmentFitTangentIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_graph(other.m_graph),
            m_report(other.m_report),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_localSize(other.m_localSize),
            m_matrix(other.m_matrix)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        SourceSegmentFitTangentIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_localSize = 0;
          m_matrix.clear();

          if (polytope.getDimension() != 1
              || polytope.getVertices().size() != 2)
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          m_localSize = polytope.getVertices().size() * sdim;
          m_matrix = computeLocalTangent(polytope, sdim);
          return *this;
        }

        Real integrate(size_t trial, size_t test) override
        {
          if (trial >= m_localSize || test >= m_localSize)
            return 0;
          return m_matrix[test * m_localSize + trial];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Faces;
        }

        SourceSegmentFitTangentIntegrator* copy() const noexcept override
        {
          return new SourceSegmentFitTangentIntegrator(*this);
        }

      private:
        std::vector<Real> computeLocalTangent(
            const Geometry::Polytope& edge,
            size_t sdim) const
        {
          const size_t localSize = edge.getVertices().size() * sdim;
          std::vector<Real> local(localSize * localSize, Real(0));
          const auto provenanceIt =
            m_report.get().interfaceEdgeProvenance.find(edge.getIndex());
          if (provenanceIt == m_report.get().interfaceEdgeProvenance.end())
            return local;

          const auto& provenance = provenanceIt->second;
          if (provenance.sourceInterfaceGraphEdge >= m_graph.get().edges.size())
            return local;

          const auto& source =
            m_graph.get().edges[provenance.sourceInterfaceGraphEdge];
          const auto& a = m_graph.get().vertices[source.v0].x;
          const auto& b = m_graph.get().vertices[source.v1].x;
          const auto ab = b - a;
          const Real denom = ab.dot(ab);
          if (denom <= Real(1e-30))
            return local;

          auto x0 = deformedVertex(m_u.get(), edge, 0, sdim);
          auto x1 = deformedVertex(m_u.get(), edge, 1, sdim);
          const Real length = (x1 - x0).norm();
          if (length <= Real(1e-14))
            return local;
          const auto edgeDirection = (Real(1) / length) * (x1 - x0);

          const auto t = (Real(1) / std::sqrt(denom)) * ab;

          static constexpr std::array<Real, 2> Points = {{
            Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
            Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
          }};
          static constexpr Real Weight = Real(0.5);

          auto projection = [&](size_t r, size_t c)
          {
            return (r == c ? Real(1) : Real(0)) - t[r] * t[c];
          };

          for (Real q : Points)
          {
            const std::array<Real, 2> N = {{Real(1) - q, q}};
            const auto x = N[0] * x0 + N[1] * x1;
            auto d = x - a;
            d -= d.dot(t) * t;
            for (size_t testNode = 0; testNode < 2; ++testNode)
            {
              for (size_t trialNode = 0; trialNode < 2; ++trialNode)
              {
                const Real lengthDerivativeSign =
                  (trialNode == 0) ? Real(-1) : Real(1);
                for (size_t r = 0; r < sdim; ++r)
                {
                  for (size_t c = 0; c < sdim; ++c)
                  {
                    const size_t row = testNode * sdim + r;
                    const size_t col = trialNode * sdim + c;
                    local[row * localSize + col] +=
                      m_weight * Weight * length
                    * N[testNode] * N[trialNode] * projection(r, c);
                    local[row * localSize + col] +=
                      m_weight * Weight
                    * N[testNode] * lengthDerivativeSign
                    * edgeDirection[c] * d[r];
                  }
                }
              }
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        std::reference_wrapper<const Geometry::InterfaceGraph> m_graph;
        std::reference_wrapper<const Geometry::LevelSetDiscretizerTrianglesReport> m_report;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localSize = 0;
        std::vector<Real> m_matrix;
    };

  template <class GridFunctionType, class ValueFunction, class GradientFunction>
  class AnalyticLevelSetFitResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LinearFormIntegratorBase<Real>;

        template <class TestFES>
        AnalyticLevelSetFitResidualIntegrator(
            const GridFunctionType& u,
            const Variational::TestFunction<TestFES>& v,
            ValueFunction value,
            GradientFunction gradient,
            Optional<Geometry::Attribute> interfaceAttribute,
            Real weight)
          : Parent(v),
            m_u(u),
            m_value(std::move(value)),
            m_gradient(std::move(gradient)),
            m_interfaceAttribute(interfaceAttribute),
            m_weight(std::max(Real(0), weight))
        {}

        AnalyticLevelSetFitResidualIntegrator(
            const AnalyticLevelSetFitResidualIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_value(other.m_value),
            m_gradient(other.m_gradient),
            m_interfaceAttribute(other.m_interfaceAttribute),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_local(other.m_local)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        AnalyticLevelSetFitResidualIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_local.clear();

          if (!isActiveInterfaceEdge(polytope))
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          m_local = computeLocalResidual(polytope, sdim);
          return *this;
        }

        Real integrate(size_t local) override
        {
          if (local >= m_local.size())
            return 0;
          return m_local[local];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Faces;
        }

        AnalyticLevelSetFitResidualIntegrator* copy() const noexcept override
        {
          return new AnalyticLevelSetFitResidualIntegrator(*this);
        }

      private:
        bool isActiveInterfaceEdge(const Geometry::Polytope& edge) const
        {
          if (edge.getDimension() != 1 || edge.getVertices().size() != 2)
            return false;
          if (!m_interfaceAttribute)
            return true;
          const auto attr = edge.getMesh().getAttribute(1, edge.getIndex());
          return attr && *attr == *m_interfaceAttribute;
        }

        std::vector<Real> computeLocalResidual(
            const Geometry::Polytope& edge,
            size_t sdim) const
        {
          const auto& fes = m_u.get().getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              edge.getDimension(), edge.getIndex());
          const size_t localSize = fe.getCount();
          std::vector<Real> local(localSize, Real(0));
          const auto coefficients =
            localEdgeCoordinateCoefficients(m_u.get(), edge);

          static constexpr std::array<Real, 2> Points = {{
            Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
            Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
          }};
          static constexpr Real Weight = Real(0.5);

          for (Real q : Points)
          {
            const Math::SpatialPoint rc{q};
            const auto x =
              deformedEdgePoint(fe, coefficients, rc, sdim);
            const auto dx =
              deformedEdgeDerivative(fe, coefficients, rc, sdim);
            const Real length = dx.norm();
            if (length <= Real(1e-14))
              continue;
            const auto direction = (Real(1) / length) * dx;
            const Real phi = m_value(x);
            const auto grad = m_gradient(x);
            for (size_t localDof = 0; localDof < localSize; ++localDof)
            {
              const size_t c = localDof % sdim;
              const Real N =
                basisComponentValue(fe, localDof, rc, c);
              const Real dN =
                basisComponentDerivative1D(fe, localDof, rc, c);
              local[localDof] += m_weight * Weight
                * (length * N * phi * grad[c]
                 + Real(0.5) * phi * phi * direction[c] * dN);
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        ValueFunction m_value;
        GradientFunction m_gradient;
        Optional<Geometry::Attribute> m_interfaceAttribute;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        std::vector<Real> m_local;
    };

  template <class GridFunctionType, class ValueFunction, class GradientFunction>
  class AnalyticLevelSetFitTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

        template <class Solution, class TrialFES, class TestFES>
        AnalyticLevelSetFitTangentIntegrator(
            const GridFunctionType& u,
            const Variational::TrialFunction<Solution, TrialFES>& du,
            const Variational::TestFunction<TestFES>& v,
            ValueFunction value,
            GradientFunction gradient,
            Optional<Geometry::Attribute> interfaceAttribute,
            Real weight)
          : Parent(du, v),
            m_u(u),
            m_value(std::move(value)),
            m_gradient(std::move(gradient)),
            m_interfaceAttribute(interfaceAttribute),
            m_weight(std::max(Real(0), weight))
        {}

        AnalyticLevelSetFitTangentIntegrator(
            const AnalyticLevelSetFitTangentIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_value(other.m_value),
            m_gradient(other.m_gradient),
            m_interfaceAttribute(other.m_interfaceAttribute),
            m_weight(other.m_weight),
            m_polytope(other.m_polytope),
            m_localSize(other.m_localSize),
            m_matrix(other.m_matrix)
        {}

        const Geometry::Polytope& getPolytope() const override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        AnalyticLevelSetFitTangentIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_localSize = 0;
          m_matrix.clear();

          if (!isActiveInterfaceEdge(polytope))
            return *this;

          const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
          const size_t sdim = mesh.getSpaceDimension();
          if (sdim != 2)
            return *this;

          const auto& fes = m_u.get().getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              polytope.getDimension(), polytope.getIndex());
          m_localSize = fe.getCount();
          m_matrix = computeLocalTangent(polytope, sdim);
          return *this;
        }

        Real integrate(size_t trial, size_t test) override
        {
          if (trial >= m_localSize || test >= m_localSize)
            return 0;
          return m_matrix[test * m_localSize + trial];
        }

        Geometry::Region getRegion() const override
        {
          return Geometry::Region::Faces;
        }

        AnalyticLevelSetFitTangentIntegrator* copy() const noexcept override
        {
          return new AnalyticLevelSetFitTangentIntegrator(*this);
        }

      private:
        bool isActiveInterfaceEdge(const Geometry::Polytope& edge) const
        {
          if (edge.getDimension() != 1 || edge.getVertices().size() != 2)
            return false;
          if (!m_interfaceAttribute)
            return true;
          const auto attr = edge.getMesh().getAttribute(1, edge.getIndex());
          return attr && *attr == *m_interfaceAttribute;
        }

        std::vector<Real> computeLocalTangent(
            const Geometry::Polytope& edge,
            size_t sdim) const
        {
          const auto& fes = m_u.get().getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              edge.getDimension(), edge.getIndex());
          const size_t localSize = fe.getCount();
          std::vector<Real> local(localSize * localSize, Real(0));
          const auto coefficients =
            localEdgeCoordinateCoefficients(m_u.get(), edge);

          static constexpr std::array<Real, 2> Points = {{
            Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
            Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
          }};
          static constexpr Real Weight = Real(0.5);

          for (Real q : Points)
          {
            const Math::SpatialPoint rc{q};
            const auto x =
              deformedEdgePoint(fe, coefficients, rc, sdim);
            const auto dx =
              deformedEdgeDerivative(fe, coefficients, rc, sdim);
            const Real length = dx.norm();
            if (length <= Real(1e-14))
              continue;
            const auto direction = (Real(1) / length) * dx;
            const Real phi = m_value(x);
            const auto grad = m_gradient(x);

            for (size_t test = 0; test < localSize; ++test)
            {
              const size_t r = test % sdim;
              const Real Ni = basisComponentValue(fe, test, rc, r);
              const Real dNi = basisComponentDerivative1D(fe, test, rc, r);
              const Real dLengthTest = direction[r] * dNi;
              const Real dPhiTest = Ni * grad[r];
              for (size_t trial = 0; trial < localSize; ++trial)
              {
                const size_t c = trial % sdim;
                const Real Nj = basisComponentValue(fe, trial, rc, c);
                const Real dNj =
                  basisComponentDerivative1D(fe, trial, rc, c);
                const Real dLengthTrial = direction[c] * dNj;
                const Real dPhiTrial = Nj * grad[c];
                const Real d2Length =
                  ((r == c ? Real(1) : Real(0)) * dNi * dNj
                 - dLengthTest * dLengthTrial) / length;

                local[test * localSize + trial] += m_weight * Weight
                  * (dLengthTrial * dPhiTest * phi
                   + length * dPhiTest * dPhiTrial
                   + phi * dPhiTrial * dLengthTest
                   + Real(0.5) * phi * phi * d2Length);
              }
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        ValueFunction m_value;
        GradientFunction m_gradient;
        Optional<Geometry::Attribute> m_interfaceAttribute;
        Real m_weight = 1;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localSize = 0;
        std::vector<Real> m_matrix;
    };

  /**
   * @brief Strict target-matrix TMOP quality term.
   *
   * The mesh coordinate map is represented by the same finite element space as
   * the displacement unknown:
   * @f[
   *   x_h = X_h + u_h, \qquad A = \nabla_{\hat x} x_h.
   * @f]
   *
   * At each element quadrature/sample point this term evaluates a fixed target
   * matrix @f$W@f$, the weighted Jacobian @f$T = A W^{-1}@f$, and the metric
   * objective
   * @f[
   *   \int_K \det(W)\,\mu(T)\,d\hat x.
   * @f]
   *
   * Current implementation status: 2D triangular vector finite element spaces
   * with two displacement components are supported and tested on P1. The public
   * API is intentionally finite-element based, so P2/H1 support can use the
   * same local-basis Jacobian path as assembly coverage grows.
   */
  template <
      class Metric = SquaredDistanceMetric,
      class TargetJacobian = IdentityTargetJacobian>
  class QualityTerm
  {
    public:
      QualityTerm(
          Metric metric = Metric{},
          TargetJacobian target = TargetJacobian{},
          Real weight = 1)
        : m_metric(std::move(metric)),
          m_target(std::move(target)),
          m_weight(std::max(Real(0), weight))
      {}

      explicit QualityTerm(Real weight)
        : QualityTerm(Metric{}, TargetJacobian{}, weight)
      {}

      QualityTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      QualityTerm& setQuadratureOrder(size_t order)
      {
        m_quadratureOrder = order;
        return *this;
      }

      Real getWeight() const
      {
        return m_weight;
      }

      size_t getQuadratureOrder() const
      {
        return m_quadratureOrder;
      }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return TMOPQualityResidualIntegrator<
          Variational::GridFunction<FES, Data>, Metric, TargetJacobian>(
              u, v, m_metric, m_target, m_weight, m_quadratureOrder);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return TMOPQualityTangentIntegrator<
          Variational::GridFunction<FES, Data>, Metric, TargetJacobian>(
              u, du, v, m_metric, m_target, m_weight, m_quadratureOrder);
      }

      template <class FES, class Data>
      Real energy(const Variational::GridFunction<FES, Data>& u) const
      {
        const auto& mesh = u.getFiniteElementSpace().getMesh();
        const auto& conn = mesh.getConnectivity();
        Real value = 0;
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          if (conn.getGeometry(2, cellIndex)
              != Geometry::Polytope::Type::Triangle)
            continue;
          auto cellIterator = mesh.getPolytope(2, cellIndex);
          const auto& cell = *cellIterator;
          const auto& fes = u.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          requireStrictTMOPCell(
              cell, fes.getVectorDimension(), fe.getCount());
          const auto coefficients = localCoordinateCoefficients(u, cell);
          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder,
              conn.getGeometry(2, cellIndex));
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const Math::SpatialMatrix<Real> A =
              deformedCoordinateJacobian(
                  fe, coefficients, rc, fes.getVectorDimension());
            const Math::SpatialMatrix<Real> W = m_target.evaluate(cell, rc);
            const Real detW = W.determinant();
            if (std::abs(detW) <= Real(1e-30))
              throw std::runtime_error(
                  "TMOP::QualityTerm target Jacobian is singular.");
            const Math::SpatialMatrix<Real> T = A * W.inverse();
            value += qf.getWeight(q) * detW * m_metric.value(T);
          }
        }
        return m_weight * value;
      }

    private:
      Metric m_metric;
      TargetJacobian m_target;
      Real m_weight = 1;
      size_t m_quadratureOrder = 2;
  };

  QualityTerm() -> QualityTerm<>;
  QualityTerm(Real) -> QualityTerm<>;
  template <class Metric, class TargetJacobian>
  QualityTerm(Metric, TargetJacobian, Real) -> QualityTerm<Metric, TargetJacobian>;
  template <class Metric, class TargetJacobian>
  QualityTerm(Metric, TargetJacobian) -> QualityTerm<Metric, TargetJacobian>;

  /**
   * @brief Legacy edge-spring smoother, not a TMOP quality metric.
   *
   * This term penalizes deviations of deformed cell edge lengths from an
   * equilateral rest length. It can remain useful as a diagnostic smoother, but
   * it is intentionally not exposed as TMOP::QualityTerm and must not be
   * reported as TMOP energy.
   */
  class EdgeSpringSmoothingTerm
  {
    public:
      explicit EdgeSpringSmoothingTerm(Real weight = 1)
        : m_weight(std::max(Real(0), weight))
      {}

      EdgeSpringSmoothingTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      Real getWeight() const
      {
        return m_weight;
      }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return EdgeSpringQualityResidualIntegrator<
          Variational::GridFunction<FES, Data>>(u, v, m_weight);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return EdgeSpringQualityTangentIntegrator<
          Variational::GridFunction<FES, Data>>(u, du, v, m_weight);
      }

    private:
      Real m_weight = 1;
  };

  /**
   * @brief Quadratic displacement-deviation term.
   *
   * Energy:
   * @f[
   *   J_{\mathrm{dev}}(u) = \frac{1}{2}\int_\Omega |u|^2\,dX.
   * @f]
   *
   * Residual:
   * @f[
   *   R(u; v) = \int_\Omega u\cdot v\,dX.
   * @f]
   *
   * Tangent:
   * @f[
   *   a(du, v) = \int_\Omega du\cdot v\,dX.
   * @f]
   */
  class DeviationTerm
  {
    public:
      explicit DeviationTerm(Real weight = 1)
        : m_weight(std::max(Real(0), weight))
      {}

      DeviationTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      Real getWeight() const
      {
        return m_weight;
      }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return DeviationResidualIntegrator<
          Variational::GridFunction<FES, Data>>(u, v, m_weight);
      }

      template <class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return DeviationTangentIntegrator(du, v, m_weight);
      }

      template <class FES, class Data>
      Real energy(const Variational::GridFunction<FES, Data>& u) const
      {
        const auto& mesh = u.getFiniteElementSpace().getMesh();
        const size_t sdim = mesh.getSpaceDimension();
        if (sdim != 2)
          return 0;

        const auto& conn = mesh.getConnectivity();
        const auto& data = u.getData();
        const Index vertexCount = static_cast<Index>(mesh.getVertexCount());
        Real value = 0;
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          if (conn.getGeometry(2, cellIndex)
              != Geometry::Polytope::Type::Triangle)
            continue;
          const auto& cell = conn.getPolytope(2, cellIndex);
          const auto x0 = mesh.getVertexCoordinates(cell(0));
          const auto x1 = mesh.getVertexCoordinates(cell(1));
          const auto x2 = mesh.getVertexCoordinates(cell(2));
          const Real area = triangleArea2D(x0, x1, x2);
          for (size_t a = 0; a < 3; ++a)
          {
            const Index va = cell(a);
            for (size_t b = 0; b < 3; ++b)
            {
              const Index vb = cell(b);
              const Real mass = area
                * ((a == b) ? Real(1) / Real(6) : Real(1) / Real(12));
              for (size_t c = 0; c < sdim; ++c)
                value += mass
                  * data(va + static_cast<Index>(c) * vertexCount)
                  * data(vb + static_cast<Index>(c) * vertexCount);
            }
          }
        }
        return Real(0.5) * m_weight * value;
      }

    private:
      Real m_weight = 1;
  };

  /**
   * @brief Source-segment fit energy for fitted level-set interfaces.
   *
   * This term evaluates and linearizes
   * @f$ \frac12\int_\Gamma \mathrm{dist}(x,\Gamma_{P1})^2\,dS @f$ on the
   * already cut interface, using LevelSetDiscretizerTriangles provenance to map
   * each output interface edge back to its source InterfaceGraph segment.
   */
  class LevelSetFitTerm
  {
    public:
      LevelSetFitTerm(
          const Geometry::InterfaceGraph& graph,
          const Geometry::LevelSetDiscretizerTrianglesReport& report,
          Real weight = 1)
        : m_graph(graph),
          m_report(report),
          m_weight(std::max(Real(0), weight))
      {}

      LevelSetFitTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      Real getWeight() const
      {
        return m_weight;
      }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return SourceSegmentFitResidualIntegrator<
          Variational::GridFunction<FES, Data>>(
              u, v, m_graph.get(), m_report.get(), m_weight);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return SourceSegmentFitTangentIntegrator<
          Variational::GridFunction<FES, Data>>(
              u, du, v, m_graph.get(), m_report.get(), m_weight);
      }

      template <class Mesh>
      Real sourceSegmentDistanceEnergy(const Mesh& fittedMesh) const
      {
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(fittedMesh, 1, 0);

        static constexpr std::array<Real, 2> Points = {{
          Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
          Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
        }};
        static constexpr Real Weight = Real(0.5);

        const auto& conn = fittedMesh.getConnectivity();
        Real res = 0;

        for (const auto& item : m_report.get().interfaceEdgeProvenance)
        {
          const Index outputEdge = item.first;
          const auto& provenance = item.second;
          if (provenance.sourceInterfaceGraphEdge >= m_graph.get().edges.size())
            continue;

          const auto& edge = conn.getPolytope(1, outputEdge);
          const auto x0 = fittedMesh.getVertexCoordinates(edge(0));
          const auto x1 = fittedMesh.getVertexCoordinates(edge(1));
          const Real length = (x1 - x0).norm();
          if (length <= Real(0))
            continue;

          const auto& source =
            m_graph.get().edges[provenance.sourceInterfaceGraphEdge];
          const auto& a = m_graph.get().vertices[source.v0].x;
          const auto& b = m_graph.get().vertices[source.v1].x;

          for (Real q : Points)
          {
            const auto x = (Real(1) - q) * x0 + q * x1;
            res += Weight * length * squaredDistanceToSegment(x, a, b);
          }
        }

        return Real(0.5) * m_weight * res;
      }

    private:
      static Real squaredDistanceToSegment(
          const Math::SpatialPoint& x,
          const Math::SpatialPoint& a,
          const Math::SpatialPoint& b)
      {
        const auto ab = b - a;
        const Real denom = ab.dot(ab);
        if (denom <= Real(0))
          return (x - a).squaredNorm();

        const Real t =
          std::clamp((x - a).dot(ab) / denom, Real(0), Real(1));
        const auto closest = a + t * ab;
        return (x - closest).squaredNorm();
      }

      std::reference_wrapper<const Geometry::InterfaceGraph> m_graph;
      std::reference_wrapper<const Geometry::LevelSetDiscretizerTrianglesReport> m_report;
      Real m_weight = 1;
  };

  /**
   * @brief Analytic level-set interface-fit penalty.
   *
   * This term evaluates the current fitted interface against an analytic level
   * set rather than against the source P1 cut segment:
   * @f[
   *   J_\phi(u) = \frac12 \int_\Gamma \phi(X + u)^2\,dS.
   * @f]
   *
   * The residual uses @f$\phi\nabla\phi@f$ on interface edges. The tangent is a
   * Gauss-Newton linearization with the interface edge-length variation
   * included, but without the analytic Hessian of @f$\phi@f$. This makes it
   * exact for affine level sets and a robust first production fit term for
   * smooth nonlinear level sets such as moving circles.
   */
  template <class ValueFunction, class GradientFunction>
  class AnalyticLevelSetFitTerm
  {
    public:
      AnalyticLevelSetFitTerm(
          ValueFunction value,
          GradientFunction gradient,
          Optional<Geometry::Attribute> interfaceAttribute = {},
          Real weight = 1)
        : m_value(std::move(value)),
          m_gradient(std::move(gradient)),
          m_interfaceAttribute(interfaceAttribute),
          m_weight(std::max(Real(0), weight))
      {}

      AnalyticLevelSetFitTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      Real getWeight() const
      {
        return m_weight;
      }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return AnalyticLevelSetFitResidualIntegrator<
          Variational::GridFunction<FES, Data>, ValueFunction, GradientFunction>(
              u, v, m_value, m_gradient, m_interfaceAttribute, m_weight);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return AnalyticLevelSetFitTangentIntegrator<
          Variational::GridFunction<FES, Data>, ValueFunction, GradientFunction>(
              u, du, v, m_value, m_gradient, m_interfaceAttribute, m_weight);
      }

      template <class FES, class Data>
      Real energy(const Variational::GridFunction<FES, Data>& u) const
      {
        const auto& mesh = u.getFiniteElementSpace().getMesh();
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        const auto& conn = mesh.getConnectivity();
        const size_t sdim = mesh.getSpaceDimension();
        if (sdim != 2)
          return 0;

        static constexpr std::array<Real, 2> Points = {{
          Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
          Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
        }};
        static constexpr Real Weight = Real(0.5);

        Real value = 0;
        for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
        {
          if (!isActiveInterfaceEdge(mesh, e))
            continue;
          const auto edgeIterator = mesh.getPolytope(1, e);
          const auto& edge = *edgeIterator;
          const auto& fes = u.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              edge.getDimension(), edge.getIndex());
          const auto coefficients = localEdgeCoordinateCoefficients(u, edge);
          for (Real q : Points)
          {
            const Math::SpatialPoint rc{q};
            const auto x = deformedEdgePoint(
                fe, coefficients, rc, fes.getVectorDimension());
            const auto dx = deformedEdgeDerivative(
                fe, coefficients, rc, fes.getVectorDimension());
            const Real length = dx.norm();
            if (length <= Real(0))
              continue;
            const Real phi = m_value(x);
            value += Weight * length * phi * phi;
          }
        }
        return Real(0.5) * m_weight * value;
      }

      template <class Mesh>
      Real energy(const Mesh& mesh) const
      {
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        static constexpr std::array<Real, 2> Points = {{
          Real(0.5) - Real(0.5) / Real(1.7320508075688772935274463415059),
          Real(0.5) + Real(0.5) / Real(1.7320508075688772935274463415059)
        }};
        static constexpr Real Weight = Real(0.5);

        const auto& conn = mesh.getConnectivity();
        Real value = 0;
        for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
        {
          if (!isActiveInterfaceEdge(mesh, e))
            continue;
          const auto edgeIterator = mesh.getPolytope(1, e);
          const auto& edge = *edgeIterator;
          const auto& trans = edge.getTransformation();
          for (Real q : Points)
          {
            const Math::SpatialPoint rc{q};
            Math::SpatialPoint x;
            trans.transform(x, rc);
            Math::SpatialMatrix<Real> J;
            trans.jacobian(J, rc);
            if (J.rows() != static_cast<std::uint8_t>(mesh.getSpaceDimension())
                || J.cols() != 1)
              continue;
            Real length2 = 0;
            for (std::uint8_t d = 0; d < J.rows(); ++d)
              length2 += J(d, 0) * J(d, 0);
            const Real length = std::sqrt(length2);
            if (length <= Real(0))
              continue;
            const Real phi = m_value(x);
            value += Weight * length * phi * phi;
          }
        }
        return Real(0.5) * m_weight * value;
      }

    private:
      template <class Mesh>
      bool isActiveInterfaceEdge(const Mesh& mesh, Index edge) const
      {
        if (!m_interfaceAttribute)
          return true;
        const auto attr = mesh.getAttribute(1, edge);
        return attr && *attr == *m_interfaceAttribute;
      }

      ValueFunction m_value;
      GradientFunction m_gradient;
      Optional<Geometry::Attribute> m_interfaceAttribute;
      Real m_weight = 1;
  };

  template <class ValueFunction, class GradientFunction>
  AnalyticLevelSetFitTerm(
      ValueFunction,
      GradientFunction,
      Optional<Geometry::Attribute>,
      Real) -> AnalyticLevelSetFitTerm<ValueFunction, GradientFunction>;

  template <class ValueFunction, class GradientFunction>
  AnalyticLevelSetFitTerm(
      ValueFunction,
      GradientFunction,
      Optional<Geometry::Attribute>) ->
        AnalyticLevelSetFitTerm<ValueFunction, GradientFunction>;

  template <class ValueFunction, class GradientFunction>
  AnalyticLevelSetFitTerm(
      ValueFunction,
      GradientFunction) ->
        AnalyticLevelSetFitTerm<ValueFunction, GradientFunction>;
}

#endif
