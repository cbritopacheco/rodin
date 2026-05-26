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

  // FES-independent: the raw displacement-field coefficients for every local
  // dof of the polytope, taken straight from the GridFunction data via the
  // FES's own local->global map. No assumption about dof layout
  // (node*vdim+component), node ordering, or basis nodality, so this is
  // correct for P1, P2 (Fekete/Dubiner) and any future element. The geometry
  // is supplied separately by the polytope transformation, never
  // reconstructed from "nodal coordinates" (which only exist for a nodal
  // interpolatory basis).
  template <class GridFunctionType>
  std::vector<Real> localDisplacementCoefficients(
      const GridFunctionType& u,
      const Geometry::Polytope& polytope)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& fe = fes.getFiniteElement(
        polytope.getDimension(), polytope.getIndex());
    const size_t localSize = fe.getCount();
    std::vector<Real> coeff(localSize, Real(0));
    const auto& data = u.getData();
    for (size_t local = 0; local < localSize; ++local)
      coeff[local] = data(fes.getGlobalIndex(
          {polytope.getDimension(), polytope.getIndex()},
          static_cast<Index>(local)));
    return coeff;
  }

  // grad_ref X_h: the polytope transformation's own reference Jacobian. Basis
  // agnostic and validated to carry P2 parametric curvature.
  inline Math::SpatialMatrix<Real> referenceGeometryJacobian(
      const Geometry::Polytope& polytope,
      const Math::SpatialPoint& rc)
  {
    return toSpatialMatrix(Geometry::Point(polytope, rc).getJacobian());
  }

  // A = grad_ref x_h = grad_ref X_h + grad_ref u_h, with grad_ref X_h from
  // the transformation (FES-independent) and grad_ref u_h = sum of the true
  // displacement dof coefficients times the FES vector-basis Jacobian. No
  // per-local nodal-coordinate reconstruction, so this is correct at any
  // order.
  template <class ElementType>
  Math::SpatialMatrix<Real> deformedCoordinateJacobian(
      const Geometry::Polytope& cell,
      const ElementType& fe,
      const std::vector<Real>& displacement,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    Math::SpatialMatrix<Real> A = referenceGeometryJacobian(cell, rc);
    for (size_t local = 0; local < displacement.size(); ++local)
      A += displacement[local] * basisJacobian2D(fe, local, rc);
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
        cell, fe, localDisplacementCoefficients(u, cell), rc,
        fes.getVectorDimension());
  }

  template <class ElementType>
  Math::SpatialPoint basisVectorValue(
      const ElementType& fe,
      size_t local,
      const Math::SpatialPoint& rc);

  template <class ElementType>
  Math::SpatialPoint deformedPoint(
      const Geometry::Polytope& polytope,
      const ElementType& fe,
      const std::vector<Real>& displacement,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    Math::SpatialPoint x = Geometry::Point(polytope, rc).getPhysicalCoordinates();
    for (size_t local = 0; local < displacement.size(); ++local)
    {
      const auto phi = basisVectorValue(fe, local, rc);
      for (size_t c = 0; c < vdim; ++c)
        x[static_cast<std::uint8_t>(c)] +=
          displacement[local] * phi[static_cast<std::uint8_t>(c)];
    }
    return x;
  }

  // The full vdim-vector value of the FES vector basis function for a local
  // dof. No component-from-layout assumption (vs local%vdim), so this is
  // correct for any vector element ordering/order.
  template <class ElementType>
  Math::SpatialPoint basisVectorValue(
      const ElementType& fe,
      size_t local,
      const Math::SpatialPoint& rc)
  {
    return fe.getBasis(local)(rc);
  }

  // d/dr of that vector basis function along the (1D) edge reference: column
  // 0 of its Jacobian.
  template <class ElementType>
  Math::SpatialPoint basisVectorDerivative1D(
      const ElementType& fe,
      size_t local,
      const Math::SpatialPoint& rc)
  {
    const auto J = toSpatialMatrix(fe.getBasis(local).getJacobian()(rc));
    Math::SpatialPoint d(J.rows());
    for (std::uint8_t c = 0; c < J.rows(); ++c)
      d[c] = J(c, 0);
    return d;
  }

  // x_h(rc) = X_h(rc) + u_h(rc): geometry from the edge transformation
  // (FES-independent, carries P2 curvature), displacement from the true dof
  // coefficients times the FES vector basis. No nodal-coordinate
  // reconstruction.
  template <class ElementType>
  Math::SpatialPoint deformedEdgePoint(
      const Geometry::Polytope& edge,
      const ElementType& fe,
      const std::vector<Real>& displacement,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    Math::SpatialPoint x = Geometry::Point(edge, rc).getPhysicalCoordinates();
    for (size_t local = 0; local < displacement.size(); ++local)
    {
      const auto phi = basisVectorValue(fe, local, rc);
      for (size_t c = 0; c < vdim; ++c)
        x[static_cast<std::uint8_t>(c)] +=
          displacement[local] * phi[static_cast<std::uint8_t>(c)];
    }
    return x;
  }

  // d x_h / dr = d X_h/dr + d u_h/dr.
  template <class ElementType>
  Math::SpatialPoint deformedEdgeDerivative(
      const Geometry::Polytope& edge,
      const ElementType& fe,
      const std::vector<Real>& displacement,
      const Math::SpatialPoint& rc,
      size_t vdim)
  {
    const auto Jx = referenceGeometryJacobian(edge, rc);
    Math::SpatialPoint dx(static_cast<std::uint8_t>(vdim));
    for (size_t c = 0; c < vdim; ++c)
      dx[static_cast<std::uint8_t>(c)] = Jx(static_cast<std::uint8_t>(c), 0);
    for (size_t local = 0; local < displacement.size(); ++local)
    {
      const auto dphi = basisVectorDerivative1D(fe, local, rc);
      for (size_t c = 0; c < vdim; ++c)
        dx[static_cast<std::uint8_t>(c)] +=
          displacement[local] * dphi[static_cast<std::uint8_t>(c)];
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
          const auto displacement =
            localDisplacementCoefficients(m_u.get(), cell);
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
            Math::SpatialMatrix<Real> A = referenceGeometryJacobian(cell, rc);
            for (size_t l = 0; l < m_localCount; ++l)
            {
              basis[l] = basisJacobian2D(fe, l, rc);
              A += displacement[l] * basis[l];
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
          const auto displacement =
            localDisplacementCoefficients(m_u.get(), cell);
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
            Math::SpatialMatrix<Real> A = referenceGeometryJacobian(cell, rc);
            for (size_t l = 0; l < m_localCount; ++l)
            {
              basis[l] = basisJacobian2D(fe, l, rc);
              A += displacement[l] * basis[l];
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
        // R_l = w * integral_K (u_h . Phi_l) dX, evaluated by quadrature on
        // the FES vector basis with the physical measure |det grad_ref X_h|.
        // FES-independent (no P1 vertex layout, no nodal coordinates), so
        // residual and tangent stay mutually consistent at any order.
        std::vector<Real> computeLocalResidual(
            const Geometry::Polytope& cell,
            size_t sdim) const
        {
          const auto& fes = m_u.get().getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          const size_t localSize = fe.getCount();
          std::vector<Real> local(localSize, Real(0));
          const auto displacement =
            localDisplacementCoefficients(m_u.get(), cell);
          const auto geometry = cell.getMesh().getGeometry(
              cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(m_quadratureOrder, geometry);
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const Real wdet = qf.getWeight(q)
              * std::abs(referenceGeometryJacobian(cell, rc).determinant());
            std::vector<Math::SpatialPoint> phi(localSize);
            Math::SpatialPoint uval(static_cast<std::uint8_t>(sdim));
            for (size_t c = 0; c < sdim; ++c)
              uval[static_cast<std::uint8_t>(c)] = 0;
            for (size_t l = 0; l < localSize; ++l)
            {
              phi[l] = basisVectorValue(fe, l, rc);
              for (size_t c = 0; c < sdim; ++c)
                uval[static_cast<std::uint8_t>(c)] +=
                  displacement[l] * phi[l][static_cast<std::uint8_t>(c)];
            }
            for (size_t l = 0; l < localSize; ++l)
            {
              Real d = 0;
              for (size_t c = 0; c < sdim; ++c)
                d += uval[static_cast<std::uint8_t>(c)]
                   * phi[l][static_cast<std::uint8_t>(c)];
              local[l] += m_weight * wdet * d;
            }
          }
          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        Real m_weight = 1;
        size_t m_quadratureOrder = 4;
        const Geometry::Polytope* m_polytope = nullptr;
        std::vector<Real> m_local;
    };

  template <class FES>
  class DeviationTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
    {
      public:
        using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

        template <class Solution, class TestFES>
        DeviationTangentIntegrator(
            const Variational::TrialFunction<Solution, FES>& du,
            const Variational::TestFunction<TestFES>& v,
            Real weight)
          : Parent(du, v),
            m_fes(&du.getFiniteElementSpace()),
            m_weight(std::max(Real(0), weight))
        {}

        DeviationTangentIntegrator(
            const DeviationTangentIntegrator& other)
          : Parent(other),
            m_fes(other.m_fes),
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

        // M_{l,m} = w * integral_K (Phi_l . Phi_m) dX, quadrature on the FES
        // vector basis with the physical measure. Same basis/quadrature as
        // the residual, so R = M u holds exactly at any element order.
        DeviationTangentIntegrator& setPolytope(
            const Geometry::Polytope& polytope) override
        {
          m_polytope = &polytope;
          m_localSize = 0;
          m_matrix.clear();

          if (polytope.getDimension() != 2
              || polytope.getVertices().size() != 3)
            return *this;

          const auto& cell = polytope;
          const auto& fe = m_fes->getFiniteElement(
              cell.getDimension(), cell.getIndex());
          const size_t sdim = m_fes->getVectorDimension();
          m_localSize = fe.getCount();
          m_matrix.assign(m_localSize * m_localSize, Real(0));

          const auto geometry = cell.getMesh().getGeometry(
              cell.getDimension(), cell.getIndex());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(m_quadratureOrder, geometry);
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const Real wdet = qf.getWeight(q)
              * std::abs(referenceGeometryJacobian(cell, rc).determinant());
            std::vector<Math::SpatialPoint> phi(m_localSize);
            for (size_t l = 0; l < m_localSize; ++l)
              phi[l] = basisVectorValue(fe, l, rc);
            for (size_t row = 0; row < m_localSize; ++row)
              for (size_t col = 0; col < m_localSize; ++col)
              {
                Real d = 0;
                for (size_t c = 0; c < sdim; ++c)
                  d += phi[row][static_cast<std::uint8_t>(c)]
                     * phi[col][static_cast<std::uint8_t>(c)];
                m_matrix[row * m_localSize + col] += m_weight * wdet * d;
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
        const FES* m_fes = nullptr;
        Real m_weight = 1;
        size_t m_quadratureOrder = 4;
        const Geometry::Polytope* m_polytope = nullptr;
        size_t m_localSize = 0;
        std::vector<Real> m_matrix;
    };

  template <class Solution, class TrialFES, class TestFES>
  DeviationTangentIntegrator(
      const Variational::TrialFunction<Solution, TrialFES>&,
      const Variational::TestFunction<TestFES>&,
      Real) -> DeviationTangentIntegrator<TrialFES>;

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
            Real weight,
            size_t quadratureOrder)
          : Parent(v),
            m_u(u),
            m_value(std::move(value)),
            m_gradient(std::move(gradient)),
            m_interfaceAttribute(interfaceAttribute),
            m_weight(std::max(Real(0), weight)),
            m_quadratureOrder(quadratureOrder)
        {}

        AnalyticLevelSetFitResidualIntegrator(
            const AnalyticLevelSetFitResidualIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_value(other.m_value),
            m_gradient(other.m_gradient),
            m_interfaceAttribute(other.m_interfaceAttribute),
            m_weight(other.m_weight),
            m_quadratureOrder(other.m_quadratureOrder),
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
          const auto displacement =
            localDisplacementCoefficients(m_u.get(), edge);

          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder, Geometry::Polytope::Type::Segment);
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const auto x =
              deformedEdgePoint(edge, fe, displacement, rc, sdim);
            const auto dx =
              deformedEdgeDerivative(edge, fe, displacement, rc, sdim);
            const Real length = dx.norm();
            if (length <= Real(1e-14))
              continue;
            const auto direction = (Real(1) / length) * dx;
            const Real phi = m_value(x);
            const auto grad = m_gradient(x);
            for (size_t localDof = 0; localDof < localSize; ++localDof)
            {
              const auto phiVec = basisVectorValue(fe, localDof, rc);
              const auto dphiVec = basisVectorDerivative1D(fe, localDof, rc);
              Real gradPhi = 0, dirDphi = 0;
              for (size_t c = 0; c < sdim; ++c)
              {
                gradPhi += grad[c] * phiVec[static_cast<std::uint8_t>(c)];
                dirDphi +=
                  direction[c] * dphiVec[static_cast<std::uint8_t>(c)];
              }
              local[localDof] += m_weight * qf.getWeight(q)
                * (length * phi * gradPhi
                 + Real(0.5) * phi * phi * dirDphi);
            }
          }

          return local;
        }

        std::reference_wrapper<const GridFunctionType> m_u;
        ValueFunction m_value;
        GradientFunction m_gradient;
        Optional<Geometry::Attribute> m_interfaceAttribute;
        Real m_weight = 1;
        size_t m_quadratureOrder = 4;
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
            Real weight,
            size_t quadratureOrder,
          std::function<Math::SpatialMatrix<Real>(const Math::SpatialPoint&)> hessian = {})
          : Parent(du, v),
            m_u(u),
            m_value(std::move(value)),
            m_gradient(std::move(gradient)),
            m_interfaceAttribute(interfaceAttribute),
            m_weight(std::max(Real(0), weight)),
            m_quadratureOrder(quadratureOrder),
            m_hessian(std::move(hessian))
        {}

        AnalyticLevelSetFitTangentIntegrator(
            const AnalyticLevelSetFitTangentIntegrator& other)
          : Parent(other),
            m_u(other.m_u),
            m_value(other.m_value),
            m_gradient(other.m_gradient),
            m_interfaceAttribute(other.m_interfaceAttribute),
            m_weight(other.m_weight),
            m_quadratureOrder(other.m_quadratureOrder),
            m_hessian(other.m_hessian),
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
          const auto displacement =
            localDisplacementCoefficients(m_u.get(), edge);

          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder, Geometry::Polytope::Type::Segment);
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const auto x =
              deformedEdgePoint(edge, fe, displacement, rc, sdim);
            const auto dx =
              deformedEdgeDerivative(edge, fe, displacement, rc, sdim);
            const Real length = dx.norm();
            if (length <= Real(1e-14))
              continue;
            const auto direction = (Real(1) / length) * dx;
            const Real phi = m_value(x);
            const auto grad = m_gradient(x);
            const bool hasHessian = static_cast<bool>(m_hessian);
            const Math::SpatialMatrix<Real> hess = hasHessian
              ? m_hessian(x)
              : Math::SpatialMatrix<Real>();

            std::vector<Math::SpatialPoint> phiVec(localSize);
            std::vector<Math::SpatialPoint> dphiVec(localSize);
            std::vector<Real> dLength(localSize, 0);
            std::vector<Real> dPhi(localSize, 0);
            for (size_t i = 0; i < localSize; ++i)
            {
              phiVec[i] = basisVectorValue(fe, i, rc);
              dphiVec[i] = basisVectorDerivative1D(fe, i, rc);
              for (size_t k = 0; k < sdim; ++k)
              {
                dLength[i] +=
                  direction[k] * dphiVec[i][static_cast<std::uint8_t>(k)];
                dPhi[i] +=
                  grad[k] * phiVec[i][static_cast<std::uint8_t>(k)];
              }
            }

            for (size_t test = 0; test < localSize; ++test)
            {
              const Real dLengthTest = dLength[test];
              const Real dPhiTest = dPhi[test];
              for (size_t trial = 0; trial < localSize; ++trial)
              {
                const Real dLengthTrial = dLength[trial];
                const Real dPhiTrial = dPhi[trial];
                Real dphiDot = 0;
                for (size_t k = 0; k < sdim; ++k)
                  dphiDot += dphiVec[test][static_cast<std::uint8_t>(k)]
                           * dphiVec[trial][static_cast<std::uint8_t>(k)];
                const Real d2Length =
                  (dphiDot - dLengthTest * dLengthTrial) / length;
                Real d2Phi = 0;
                if (hasHessian)
                {
                  for (size_t a = 0; a < sdim; ++a)
                    for (size_t b = 0; b < sdim; ++b)
                      d2Phi += phiVec[test][static_cast<std::uint8_t>(a)]
                            * hess(a, b)
                            * phiVec[trial][static_cast<std::uint8_t>(b)];
                }

                local[test * localSize + trial] += m_weight * qf.getWeight(q)
                  * (dLengthTrial * dPhiTest * phi
                   + length * dPhiTest * dPhiTrial
                   + length * phi * d2Phi
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
        size_t m_quadratureOrder = 4;
        std::function<Math::SpatialMatrix<Real>(const Math::SpatialPoint&)> m_hessian;
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
          const auto displacement = localDisplacementCoefficients(u, cell);
          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder,
              conn.getGeometry(2, cellIndex));
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const Math::SpatialMatrix<Real> A =
              deformedCoordinateJacobian(
                  cell, fe, displacement, rc, fes.getVectorDimension());
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
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t sdim = mesh.getSpaceDimension();
        if (sdim != 2)
          return 0;

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
          const auto& fe = fes.getFiniteElement(
              cell.getDimension(), cell.getIndex());
          const size_t localSize = fe.getCount();
          const auto displacement = localDisplacementCoefficients(u, cell);
          const auto& qf = QF::PolytopeQuadratureFormula::get(
              4, conn.getGeometry(2, cellIndex));
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const Real wdet = qf.getWeight(q)
              * std::abs(referenceGeometryJacobian(cell, rc).determinant());
            Math::SpatialPoint uval(static_cast<std::uint8_t>(sdim));
            for (size_t c = 0; c < sdim; ++c)
              uval[static_cast<std::uint8_t>(c)] = 0;
            for (size_t l = 0; l < localSize; ++l)
            {
              const auto phi = basisVectorValue(fe, l, rc);
              for (size_t c = 0; c < sdim; ++c)
                uval[static_cast<std::uint8_t>(c)] +=
                  displacement[l] * phi[static_cast<std::uint8_t>(c)];
            }
            Real sq = 0;
            for (size_t c = 0; c < sdim; ++c)
              sq += uval[static_cast<std::uint8_t>(c)]
                  * uval[static_cast<std::uint8_t>(c)];
            value += wdet * sq;
          }
        }
        return Real(0.5) * m_weight * value;
      }

    private:
      Real m_weight = 1;
  };

  template <class GridFunctionType, class ValueFunction, class GradientFunction>
  class VolumetricPhaseConsistencyResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using Parent = Variational::LinearFormIntegratorBase<Real>;

      template <class TestFES>
      VolumetricPhaseConsistencyResidualIntegrator(
          const GridFunctionType& u,
          const Variational::TestFunction<TestFES>& v,
          ValueFunction value,
          GradientFunction gradient,
          Geometry::Attribute negativeAttribute,
          Geometry::Attribute positiveAttribute,
          Real weight,
          Real epsilon,
          Real margin,
          Real normalization,
          size_t quadratureOrder)
        : Parent(v),
          m_u(u),
          m_value(std::move(value)),
          m_gradient(std::move(gradient)),
          m_negativeAttribute(negativeAttribute),
          m_positiveAttribute(positiveAttribute),
          m_weight(std::max(Real(0), weight)),
          m_epsilon(std::max(epsilon, Real(1e-12))),
          m_margin(margin),
          m_normalization(std::max(normalization, Real(1e-12))),
          m_quadratureOrder(quadratureOrder)
      {}

      VolumetricPhaseConsistencyResidualIntegrator(
          const VolumetricPhaseConsistencyResidualIntegrator& other)
        : Parent(other),
          m_u(other.m_u),
          m_value(other.m_value),
          m_gradient(other.m_gradient),
          m_negativeAttribute(other.m_negativeAttribute),
          m_positiveAttribute(other.m_positiveAttribute),
          m_weight(other.m_weight),
          m_epsilon(other.m_epsilon),
          m_margin(other.m_margin),
          m_normalization(other.m_normalization),
          m_quadratureOrder(other.m_quadratureOrder),
          m_polytope(other.m_polytope),
          m_localCount(other.m_localCount),
          m_scale(other.m_scale),
          m_basisValue(other.m_basisValue),
          m_gradPhi(other.m_gradPhi)
      {}

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      VolumetricPhaseConsistencyResidualIntegrator& setPolytope(
          const Geometry::Polytope& polytope) override
      {
        m_polytope = &polytope;
        const auto& cell = polytope;
        const auto attr = cell.getMesh().getAttribute(
            cell.getDimension(), cell.getIndex());
        const Real phaseSign = getPhaseSign(attr);
        m_localCount = 0;
        m_scale.clear();
        m_basisValue.clear();
        m_gradPhi.clear();
        if (phaseSign == Real(0))
          return *this;

        const auto& fes = m_u.get().getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(
            cell.getDimension(), cell.getIndex());
        m_localCount = fe.getCount();
        const auto displacement = localDisplacementCoefficients(m_u.get(), cell);
        const auto geometry =
          cell.getMesh().getGeometry(cell.getDimension(), cell.getIndex());
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(m_quadratureOrder, geometry);
        const size_t nq = qf.getSize();
        m_scale.assign(nq, 0);
        m_basisValue.assign(
            nq, std::vector<Math::SpatialPoint>(m_localCount));
        m_gradPhi.assign(
            nq, Math::SpatialPoint(
                static_cast<std::uint8_t>(fes.getVectorDimension())));
        for (size_t q = 0; q < nq; ++q)
        {
          const auto& rc = qf.getPoint(q);
          Math::SpatialMatrix<Real> J0;
          cell.getTransformation().jacobian(J0, rc);
          const Real wdet = qf.getWeight(q) * std::abs(J0.determinant());
          const auto x = deformedPoint(
              cell, fe, displacement, rc, fes.getVectorDimension());
          const Real z = phaseSign * m_value(x) / m_epsilon;
          if (z >= m_margin)
            continue;
          const Real psiPrime = z - m_margin;
          m_scale[q] = wdet * psiPrime * phaseSign / m_epsilon;
          m_gradPhi[q] = m_gradient(x);
          for (size_t l = 0; l < m_localCount; ++l)
            m_basisValue[q][l] = basisVectorValue(fe, l, rc);
        }
        return *this;
      }

      Real integrate(size_t local) override
      {
        if (local >= m_localCount)
          return 0;
        Real value = 0;
        for (size_t q = 0; q < m_scale.size(); ++q)
        {
          if (m_scale[q] == Real(0))
            continue;
          value += m_scale[q] * Math::dot(m_gradPhi[q], m_basisValue[q][local]);
        }
        return (m_weight / m_normalization) * value;
      }

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      VolumetricPhaseConsistencyResidualIntegrator* copy() const noexcept override
      {
        return new VolumetricPhaseConsistencyResidualIntegrator(*this);
      }

    private:
      Real getPhaseSign(const Optional<Geometry::Attribute>& attr) const
      {
        if (!attr)
          return Real(0);
        if (*attr == m_negativeAttribute)
          return Real(-1);
        if (*attr == m_positiveAttribute)
          return Real(1);
        return Real(0);
      }

      std::reference_wrapper<const GridFunctionType> m_u;
      ValueFunction m_value;
      GradientFunction m_gradient;
      Geometry::Attribute m_negativeAttribute;
      Geometry::Attribute m_positiveAttribute;
      Real m_weight = 1;
      Real m_epsilon = 1;
      Real m_margin = 1;
      Real m_normalization = 1;
      size_t m_quadratureOrder = 2;
      const Geometry::Polytope* m_polytope = nullptr;
      size_t m_localCount = 0;
      std::vector<Real> m_scale;
      std::vector<std::vector<Math::SpatialPoint>> m_basisValue;
      std::vector<Math::SpatialPoint> m_gradPhi;
  };

  template <class GridFunctionType, class ValueFunction, class GradientFunction>
  class VolumetricPhaseConsistencyTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

      template <class Solution, class TrialFES, class TestFES>
      VolumetricPhaseConsistencyTangentIntegrator(
          const GridFunctionType& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v,
          ValueFunction value,
          GradientFunction gradient,
          Geometry::Attribute negativeAttribute,
          Geometry::Attribute positiveAttribute,
          Real weight,
          Real epsilon,
          Real margin,
          Real normalization,
          size_t quadratureOrder,
          std::function<Math::SpatialMatrix<Real>(const Math::SpatialPoint&)> hessian = {})
        : Parent(du, v),
          m_u(u),
          m_value(std::move(value)),
          m_gradient(std::move(gradient)),
          m_negativeAttribute(negativeAttribute),
          m_positiveAttribute(positiveAttribute),
          m_weight(std::max(Real(0), weight)),
          m_epsilon(std::max(epsilon, Real(1e-12))),
          m_margin(margin),
          m_normalization(std::max(normalization, Real(1e-12))),
          m_quadratureOrder(quadratureOrder),
          m_hessian(std::move(hessian))
      {}

      VolumetricPhaseConsistencyTangentIntegrator(
          const VolumetricPhaseConsistencyTangentIntegrator& other)
        : Parent(other),
          m_u(other.m_u),
          m_value(other.m_value),
          m_gradient(other.m_gradient),
          m_negativeAttribute(other.m_negativeAttribute),
          m_positiveAttribute(other.m_positiveAttribute),
          m_weight(other.m_weight),
          m_epsilon(other.m_epsilon),
          m_margin(other.m_margin),
          m_normalization(other.m_normalization),
          m_quadratureOrder(other.m_quadratureOrder),
          m_hessian(other.m_hessian),
          m_polytope(other.m_polytope),
          m_localCount(other.m_localCount),
          m_scale(other.m_scale),
          m_curvatureScale(other.m_curvatureScale),
          m_basisValue(other.m_basisValue),
          m_gradPhi(other.m_gradPhi),
          m_hessianValue(other.m_hessianValue)
      {}

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      VolumetricPhaseConsistencyTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) override
      {
        m_polytope = &polytope;
        const auto& cell = polytope;
        const auto attr = cell.getMesh().getAttribute(
            cell.getDimension(), cell.getIndex());
        const Real phaseSign = getPhaseSign(attr);
        m_localCount = 0;
        m_scale.clear();
        m_curvatureScale.clear();
        m_basisValue.clear();
        m_gradPhi.clear();
        m_hessianValue.clear();
        if (phaseSign == Real(0))
          return *this;

        const auto& fes = m_u.get().getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(
            cell.getDimension(), cell.getIndex());
        m_localCount = fe.getCount();
        const auto displacement = localDisplacementCoefficients(m_u.get(), cell);
        const auto geometry =
          cell.getMesh().getGeometry(cell.getDimension(), cell.getIndex());
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(m_quadratureOrder, geometry);
        const size_t nq = qf.getSize();
        m_scale.assign(nq, 0);
        m_curvatureScale.assign(nq, 0);
        m_basisValue.assign(
            nq, std::vector<Math::SpatialPoint>(m_localCount));
        m_gradPhi.assign(
            nq, Math::SpatialPoint(
                static_cast<std::uint8_t>(fes.getVectorDimension())));
        m_hessianValue.assign(nq, Math::SpatialMatrix<Real>());
        for (size_t q = 0; q < nq; ++q)
        {
          const auto& rc = qf.getPoint(q);
          Math::SpatialMatrix<Real> J0;
          cell.getTransformation().jacobian(J0, rc);
          const Real wdet = qf.getWeight(q) * std::abs(J0.determinant());
          const auto x = deformedPoint(
              cell, fe, displacement, rc, fes.getVectorDimension());
          const Real z = phaseSign * m_value(x) / m_epsilon;
          if (z >= m_margin)
            continue;
          const Real psiPrime = z - m_margin;
          const Real psiSecond = Real(1);
          const Real factor = phaseSign / m_epsilon;
          m_scale[q] = wdet * psiSecond * factor * factor;
          m_gradPhi[q] = m_gradient(x);
          if (m_hessian)
          {
            m_curvatureScale[q] = wdet * psiPrime * factor;
            m_hessianValue[q] = m_hessian(x);
          }
          for (size_t l = 0; l < m_localCount; ++l)
            m_basisValue[q][l] = basisVectorValue(fe, l, rc);
        }
        return *this;
      }

      Real integrate(size_t trial, size_t test) override
      {
        if (trial >= m_localCount || test >= m_localCount)
          return 0;
        Real value = 0;
        for (size_t q = 0; q < m_scale.size(); ++q)
        {
          if (m_scale[q] == Real(0))
            continue;
          const Real dTrial =
            Math::dot(m_gradPhi[q], m_basisValue[q][trial]);
          const Real dTest =
            Math::dot(m_gradPhi[q], m_basisValue[q][test]);
          value += m_scale[q] * dTrial * dTest;
          if (m_curvatureScale[q] != Real(0))
          {
            Real curvature = 0;
            const auto& h = m_hessianValue[q];
            for (std::uint8_t a = 0; a < h.rows(); ++a)
              for (std::uint8_t b = 0; b < h.cols(); ++b)
                curvature += m_basisValue[q][test][a]
                           * h(a, b)
                           * m_basisValue[q][trial][b];
            value += m_curvatureScale[q] * curvature;
          }
        }
        return (m_weight / m_normalization) * value;
      }

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      VolumetricPhaseConsistencyTangentIntegrator* copy() const noexcept override
      {
        return new VolumetricPhaseConsistencyTangentIntegrator(*this);
      }

    private:
      Real getPhaseSign(const Optional<Geometry::Attribute>& attr) const
      {
        if (!attr)
          return Real(0);
        if (*attr == m_negativeAttribute)
          return Real(-1);
        if (*attr == m_positiveAttribute)
          return Real(1);
        return Real(0);
      }

      std::reference_wrapper<const GridFunctionType> m_u;
      ValueFunction m_value;
      GradientFunction m_gradient;
      Geometry::Attribute m_negativeAttribute;
      Geometry::Attribute m_positiveAttribute;
      Real m_weight = 1;
      Real m_epsilon = 1;
      Real m_margin = 1;
      Real m_normalization = 1;
      size_t m_quadratureOrder = 2;
      std::function<Math::SpatialMatrix<Real>(const Math::SpatialPoint&)> m_hessian;
      const Geometry::Polytope* m_polytope = nullptr;
      size_t m_localCount = 0;
      std::vector<Real> m_scale;
      std::vector<Real> m_curvatureScale;
      std::vector<std::vector<Math::SpatialPoint>> m_basisValue;
      std::vector<Math::SpatialPoint> m_gradPhi;
      std::vector<Math::SpatialMatrix<Real>> m_hessianValue;
  };

  template <class ValueFunction, class GradientFunction>
  class VolumetricPhaseConsistencyTerm
  {
    public:
      using HessianFunction = std::function<Math::SpatialMatrix<Real>(const Math::SpatialPoint&)>;

      VolumetricPhaseConsistencyTerm(
          ValueFunction value,
          GradientFunction gradient,
          Geometry::Attribute negativeAttribute,
          Geometry::Attribute positiveAttribute,
          Real weight = 1)
        : m_value(std::move(value)),
          m_gradient(std::move(gradient)),
          m_negativeAttribute(negativeAttribute),
          m_positiveAttribute(positiveAttribute),
          m_weight(std::max(Real(0), weight))
      {}

      template <class HessianLike>
      VolumetricPhaseConsistencyTerm(
          ValueFunction value,
          GradientFunction gradient,
          HessianLike hessian,
          Geometry::Attribute negativeAttribute,
          Geometry::Attribute positiveAttribute,
          Real weight = 1)
        : m_value(std::move(value)),
          m_gradient(std::move(gradient)),
          m_hessian(HessianFunction(std::move(hessian))),
          m_negativeAttribute(negativeAttribute),
          m_positiveAttribute(positiveAttribute),
          m_weight(std::max(Real(0), weight))
      {}

      VolumetricPhaseConsistencyTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      VolumetricPhaseConsistencyTerm& setEpsilon(Real epsilon)
      {
        m_epsilon = std::max(epsilon, Real(1e-12));
        return *this;
      }

      VolumetricPhaseConsistencyTerm& setMargin(Real margin)
      {
        m_margin = margin;
        return *this;
      }

      VolumetricPhaseConsistencyTerm& setNormalization(Real normalization)
      {
        m_normalization = normalization > Real(0) ? normalization : Real(1);
        return *this;
      }

      VolumetricPhaseConsistencyTerm& setQuadratureOrder(size_t quadratureOrder)
      {
        m_quadratureOrder = quadratureOrder;
        return *this;
      }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return VolumetricPhaseConsistencyResidualIntegrator<
          Variational::GridFunction<FES, Data>, ValueFunction, GradientFunction>(
              u, v, m_value, m_gradient,
              m_negativeAttribute, m_positiveAttribute,
              m_weight, m_epsilon, m_margin,
              m_normalization, m_quadratureOrder);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return VolumetricPhaseConsistencyTangentIntegrator<
          Variational::GridFunction<FES, Data>, ValueFunction, GradientFunction>(
              u, du, v, m_value, m_gradient,
              m_negativeAttribute, m_positiveAttribute,
              m_weight, m_epsilon, m_margin,
              m_normalization, m_quadratureOrder, m_hessian);
      }

      template <class FES, class Data>
      Real energy(const Variational::GridFunction<FES, Data>& u) const
      {
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& conn = mesh.getConnectivity();
        Real value = 0;
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          const auto attr = mesh.getAttribute(2, cellIndex);
          const Real phaseSign = getPhaseSign(attr);
          if (phaseSign == Real(0))
            continue;
          const auto cell = mesh.getPolytope(2, cellIndex);
          const auto& fe = fes.getFiniteElement(
              cell->getDimension(), cell->getIndex());
          const auto displacement = localDisplacementCoefficients(u, *cell);
          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder, conn.getGeometry(2, cellIndex));
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            Math::SpatialMatrix<Real> J0;
            cell->getTransformation().jacobian(J0, rc);
            const Real wdet = qf.getWeight(q) * std::abs(J0.determinant());
            const auto x = deformedPoint(
                *cell, fe, displacement, rc, fes.getVectorDimension());
            const Real z = phaseSign * m_value(x) / m_epsilon;
            const Real gap = std::max(Real(0), m_margin - z);
            value += Real(0.5) * wdet * gap * gap;
          }
        }
        return (m_weight / m_normalization) * value;
      }

      template <class Mesh>
      Real energy(const Mesh& mesh) const
      {
        const auto& conn = mesh.getConnectivity();
        Real value = 0;
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          const auto attr = mesh.getAttribute(2, cellIndex);
          const Real phaseSign = getPhaseSign(attr);
          if (phaseSign == Real(0))
            continue;
          const auto cell = mesh.getPolytope(2, cellIndex);
          const auto& trans = cell->getTransformation();
          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder, conn.getGeometry(2, cellIndex));
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            Math::SpatialPoint x;
            trans.transform(x, qf.getPoint(q));
            Math::SpatialMatrix<Real> J;
            trans.jacobian(J, qf.getPoint(q));
            if (J.rows() != 2 || J.cols() != 2)
              continue;
            const Real wdet = qf.getWeight(q) * std::abs(J.determinant());
            const Real z = phaseSign * m_value(x) / m_epsilon;
            const Real gap = std::max(Real(0), m_margin - z);
            value += Real(0.5) * wdet * gap * gap;
          }
        }
        return (m_weight / m_normalization) * value;
      }

      template <class Mesh>
      Index countWrongSideQuadrature(const Mesh& mesh) const
      {
        Index count = 0;
        const auto& conn = mesh.getConnectivity();
        for (Index cellIndex = 0;
             cellIndex < static_cast<Index>(mesh.getCellCount());
             ++cellIndex)
        {
          const auto attr = mesh.getAttribute(2, cellIndex);
          const Real phaseSign = getPhaseSign(attr);
          if (phaseSign == Real(0))
            continue;
          const auto cell = mesh.getPolytope(2, cellIndex);
          const auto& trans = cell->getTransformation();
          const auto& qf = QF::PolytopeQuadratureFormula::get(
              m_quadratureOrder, conn.getGeometry(2, cellIndex));
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            Math::SpatialPoint x;
            trans.transform(x, qf.getPoint(q));
            const Real z = phaseSign * m_value(x) / m_epsilon;
            if (z < Real(0))
              ++count;
          }
        }
        return count;
      }

    private:
      Real getPhaseSign(const Optional<Geometry::Attribute>& attr) const
      {
        if (!attr)
          return Real(0);
        if (*attr == m_negativeAttribute)
          return Real(-1);
        if (*attr == m_positiveAttribute)
          return Real(1);
        return Real(0);
      }

      ValueFunction m_value;
      GradientFunction m_gradient;
      HessianFunction m_hessian;
      Geometry::Attribute m_negativeAttribute;
      Geometry::Attribute m_positiveAttribute;
      Real m_weight = 1;
      Real m_epsilon = 1;
      Real m_margin = 1;
      Real m_normalization = 1;
      size_t m_quadratureOrder = 4;
  };

  template <class ValueFunction, class GradientFunction, class HessianFunction>
  VolumetricPhaseConsistencyTerm(
      ValueFunction,
      GradientFunction,
      HessianFunction,
      Geometry::Attribute,
      Geometry::Attribute,
      Real) -> VolumetricPhaseConsistencyTerm<ValueFunction, GradientFunction>;

  /**
   * @brief Analytic level-set interface-fit penalty.
   *
   * This term evaluates the current fitted interface against an analytic level
   * set rather than against the source P1 cut segment:
   * @f[
   *   J_\phi(u) = \frac12 \int_\Gamma \phi(X + u)^2\,dS.
   * @f]
   *
  * The residual uses @f$\phi\nabla\phi@f$ on interface edges. The tangent is
  * exact when an analytic Hessian is supplied, and otherwise falls back to a
  * Gauss-Newton linearization without the level-set Hessian. This keeps the
  * term usable for smooth nonlinear level sets while still supporting exact
  * Newton behavior for problems such as the moving-circle example.
   */
  template <class ValueFunction, class GradientFunction>
  class AnalyticLevelSetFitTerm
  {
    public:
      using HessianFunction = std::function<Math::SpatialMatrix<Real>(const Math::SpatialPoint&)>;

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

      AnalyticLevelSetFitTerm(
          ValueFunction value,
          GradientFunction gradient,
          HessianFunction hessian,
          Optional<Geometry::Attribute> interfaceAttribute,
          Real weight = 1)
        : m_value(std::move(value)),
          m_gradient(std::move(gradient)),
          m_hessian(std::move(hessian)),
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

      /// Normalization constant c_sigma (Knupp-Kolev-Mittal-Tomov 2021): the
      /// penalty acts with effective weight weight / c_sigma. Set c_sigma to
      /// the interface measure (total interface length in 2D) so the fit
      /// weight is portable across mesh resolutions -- the penalty then
      /// measures the *mean-square* level-set violation along the interface
      /// rather than its (resolution-dependent) integral. Default 1 (off).
      AnalyticLevelSetFitTerm& setNormalization(Real cSigma)
      {
        m_normalization = (cSigma > Real(0)) ? cSigma : Real(1);
        return *this;
      }

      Real getNormalization() const
      {
        return m_normalization;
      }

      AnalyticLevelSetFitTerm& setQuadratureOrder(size_t quadratureOrder)
      {
        m_quadratureOrder = quadratureOrder;
        return *this;
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
        return AnalyticLevelSetFitResidualIntegrator<
              Variational::GridFunction<FES, Data>, ValueFunction, GradientFunction>(
              u, v, m_value, m_gradient, m_interfaceAttribute,
              m_weight / m_normalization, m_quadratureOrder);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return AnalyticLevelSetFitTangentIntegrator<
              Variational::GridFunction<FES, Data>, ValueFunction, GradientFunction>(
              u, du, v, m_value, m_gradient, m_interfaceAttribute,
              m_weight / m_normalization, m_quadratureOrder, m_hessian);
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

        const auto& qf = QF::PolytopeQuadratureFormula::get(
            m_quadratureOrder, Geometry::Polytope::Type::Segment);

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
          const auto displacement = localDisplacementCoefficients(u, edge);
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
            const auto x = deformedEdgePoint(
                edge, fe, displacement, rc, fes.getVectorDimension());
            const auto dx = deformedEdgeDerivative(
                edge, fe, displacement, rc, fes.getVectorDimension());
            const Real length = dx.norm();
            if (length <= Real(0))
              continue;
            const Real phi = m_value(x);
            value += qf.getWeight(q) * length * phi * phi;
          }
        }
        return Real(0.5) * (m_weight / m_normalization) * value;
      }

      template <class Mesh>
      Real energy(const Mesh& mesh) const
      {
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        const auto& conn = mesh.getConnectivity();
        const auto& qf = QF::PolytopeQuadratureFormula::get(
            m_quadratureOrder, Geometry::Polytope::Type::Segment);
        Real value = 0;
        for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
        {
          if (!isActiveInterfaceEdge(mesh, e))
            continue;
          const auto edgeIterator = mesh.getPolytope(1, e);
          const auto& edge = *edgeIterator;
          const auto& trans = edge.getTransformation();
          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc = qf.getPoint(q);
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
            value += qf.getWeight(q) * length * phi * phi;
          }
        }
        return Real(0.5) * (m_weight / m_normalization) * value;
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
        HessianFunction m_hessian;
      Optional<Geometry::Attribute> m_interfaceAttribute;
      Real m_weight = 1;
      Real m_normalization = 1;
      size_t m_quadratureOrder = 4;
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

  // ── InterfaceNormalAdvectionFitTerm ──────────────────────────────────────
  //
  // Penalty: rho/2 * integral_{Gamma_h} (u·n − target(x))^2 ds
  //
  // n = 90° CCW rotation of the undeformed edge tangent.
  // target(x) = Δt * V_n(x) where V_n is the signed normal velocity.
  // The term is quadratic in u; geometry is evaluated at the undeformed
  // configuration so the tangent matrix is symmetric and constant per step.

  template <class GridFunctionType, class TargetFunction>
  class InterfaceNormalAdvectionFitResidualIntegrator final
      : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using Parent = Variational::LinearFormIntegratorBase<Real>;

      template <class TestFES>
      InterfaceNormalAdvectionFitResidualIntegrator(
          const GridFunctionType& u,
          const Variational::TestFunction<TestFES>& v,
          TargetFunction target,
          Optional<Geometry::Attribute> interfaceAttribute,
          Real weight,
          size_t quadratureOrder)
        : Parent(v),
          m_u(u),
          m_target(std::move(target)),
          m_interfaceAttribute(interfaceAttribute),
          m_weight(std::max(Real(0), weight)),
          m_quadratureOrder(quadratureOrder)
      {}

      InterfaceNormalAdvectionFitResidualIntegrator(
          const InterfaceNormalAdvectionFitResidualIntegrator& other)
        : Parent(other),
          m_u(other.m_u),
          m_target(other.m_target),
          m_interfaceAttribute(other.m_interfaceAttribute),
          m_weight(other.m_weight),
          m_quadratureOrder(other.m_quadratureOrder),
          m_polytope(other.m_polytope),
          m_local(other.m_local)
      {}

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      InterfaceNormalAdvectionFitResidualIntegrator& setPolytope(
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

      InterfaceNormalAdvectionFitResidualIntegrator* copy() const noexcept override
      {
        return new InterfaceNormalAdvectionFitResidualIntegrator(*this);
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
        const auto& fe  = fes.getFiniteElement(edge.getDimension(), edge.getIndex());
        const size_t localSize = fe.getCount();
        std::vector<Real> local(localSize, Real(0));
        const auto displacement = localDisplacementCoefficients(m_u.get(), edge);

        const auto& qf = QF::PolytopeQuadratureFormula::get(
            m_quadratureOrder, Geometry::Polytope::Type::Segment);

        for (size_t q = 0; q < qf.getSize(); ++q)
        {
          const auto& rc    = qf.getPoint(q);
          const auto   J    = referenceGeometryJacobian(edge, rc);
          const Real arcLen = std::hypot(J(0, 0), J(1, 0));
          if (arcLen <= Real(1e-14))
            continue;

          // Outward normal: 90° CCW rotation of tangent (t -> n = (-ty, tx)).
          // Sign convention: orientation-consistent with the edge parametrization.
          const Math::SpatialPoint n{ -J(1, 0) / arcLen, J(0, 0) / arcLen };

          Math::SpatialPoint x;
          edge.getTransformation().transform(x, rc);
          const Real tgt = m_target(x);

          // u_h · n at this QP
          Real uDotN = Real(0);
          for (size_t i = 0; i < localSize; ++i)
          {
            const auto phiI = basisVectorValue(fe, i, rc);
            Real p = Real(0);
            for (size_t c = 0; c < sdim; ++c)
              p += phiI[static_cast<std::uint8_t>(c)] * n[c];
            uDotN += displacement[i] * p;
          }

          for (size_t j = 0; j < localSize; ++j)
          {
            const auto phiJ = basisVectorValue(fe, j, rc);
            Real p = Real(0);
            for (size_t c = 0; c < sdim; ++c)
              p += phiJ[static_cast<std::uint8_t>(c)] * n[c];
            local[j] += m_weight * qf.getWeight(q) * arcLen * (uDotN - tgt) * p;
          }
        }
        return local;
      }

      std::reference_wrapper<const GridFunctionType> m_u;
      TargetFunction m_target;
      Optional<Geometry::Attribute> m_interfaceAttribute;
      Real   m_weight         = 1;
      size_t m_quadratureOrder = 4;
      const Geometry::Polytope* m_polytope = nullptr;
      std::vector<Real> m_local;
  };

  template <class GridFunctionType, class TargetFunction>
  class InterfaceNormalAdvectionFitTangentIntegrator final
      : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

      template <class Solution, class TrialFES, class TestFES>
      InterfaceNormalAdvectionFitTangentIntegrator(
          const GridFunctionType& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v,
          Optional<Geometry::Attribute> interfaceAttribute,
          Real weight,
          size_t quadratureOrder)
        : Parent(du, v),
          m_u(u),
          m_interfaceAttribute(interfaceAttribute),
          m_weight(std::max(Real(0), weight)),
          m_quadratureOrder(quadratureOrder)
      {}

      InterfaceNormalAdvectionFitTangentIntegrator(
          const InterfaceNormalAdvectionFitTangentIntegrator& other)
        : Parent(other),
          m_u(other.m_u),
          m_interfaceAttribute(other.m_interfaceAttribute),
          m_weight(other.m_weight),
          m_quadratureOrder(other.m_quadratureOrder),
          m_polytope(other.m_polytope),
          m_localSize(other.m_localSize),
          m_matrix(other.m_matrix)
      {}

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      InterfaceNormalAdvectionFitTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) override
      {
        m_polytope  = &polytope;
        m_localSize = 0;
        m_matrix.clear();
        if (!isActiveInterfaceEdge(polytope))
          return *this;
        const auto& mesh = m_u.get().getFiniteElementSpace().getMesh();
        const size_t sdim = mesh.getSpaceDimension();
        if (sdim != 2)
          return *this;
        const auto& fes = m_u.get().getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(polytope.getDimension(), polytope.getIndex());
        m_localSize = fe.getCount();
        m_matrix    = computeLocalTangent(polytope, sdim);
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

      InterfaceNormalAdvectionFitTangentIntegrator* copy() const noexcept override
      {
        return new InterfaceNormalAdvectionFitTangentIntegrator(*this);
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
        const auto& fe  = fes.getFiniteElement(edge.getDimension(), edge.getIndex());
        const size_t localSize = fe.getCount();
        std::vector<Real> local(localSize * localSize, Real(0));

        const auto& qf = QF::PolytopeQuadratureFormula::get(
            m_quadratureOrder, Geometry::Polytope::Type::Segment);

        for (size_t q = 0; q < qf.getSize(); ++q)
        {
          const auto& rc    = qf.getPoint(q);
          const auto   J    = referenceGeometryJacobian(edge, rc);
          const Real arcLen = std::hypot(J(0, 0), J(1, 0));
          if (arcLen <= Real(1e-14))
            continue;
          const Math::SpatialPoint n{ -J(1, 0) / arcLen, J(0, 0) / arcLen };

          // phi_i · n for all local DOFs
          std::vector<Real> phiDotN(localSize, Real(0));
          for (size_t i = 0; i < localSize; ++i)
          {
            const auto phi = basisVectorValue(fe, i, rc);
            for (size_t c = 0; c < sdim; ++c)
              phiDotN[i] += phi[static_cast<std::uint8_t>(c)] * n[c];
          }

          for (size_t test = 0; test < localSize; ++test)
            for (size_t trial = 0; trial < localSize; ++trial)
              local[test * localSize + trial] +=
                m_weight * qf.getWeight(q) * arcLen * phiDotN[test] * phiDotN[trial];
        }
        return local;
      }

      std::reference_wrapper<const GridFunctionType> m_u;
      Optional<Geometry::Attribute> m_interfaceAttribute;
      Real   m_weight          = 1;
      size_t m_quadratureOrder = 4;
      const Geometry::Polytope* m_polytope = nullptr;
      size_t m_localSize = 0;
      std::vector<Real> m_matrix;
  };

  /**
   * @brief Advection-driven interface fit term.
   *
   * Assembles the quadratic penalty:
   * @f[
   *   \frac{\rho}{2} \int_{\Gamma_h}
   *     \bigl(u \cdot n - \Delta t\,V_n(x)\bigr)^2 \,\mathrm{d}s
   * @f]
   * where @f$n@f$ is the 90° CCW rotation of the undeformed edge tangent and
   * @f$V_n(x) = V(x)\cdot n(x)@f$ is the signed normal velocity supplied
   * by the caller via @p target(x) = Δt * V_n(x).
   *
   * Since @f$n@f$ is evaluated at the current (undeformed) geometry, the term
   * is linear in @f$u@f$ (residual) and has a symmetric positive semi-definite
   * tangent that is constant per time step — cheaper to assemble than the
   * nonlinear @f$\phi(x+u)@f$ term.
   */
  template <class TargetFunction>
  class InterfaceNormalAdvectionFitTerm
  {
    public:
      InterfaceNormalAdvectionFitTerm(
          TargetFunction target,
          Optional<Geometry::Attribute> interfaceAttribute = {},
          Real weight = 1)
        : m_target(std::move(target)),
          m_interfaceAttribute(interfaceAttribute),
          m_weight(std::max(Real(0), weight))
      {}

      InterfaceNormalAdvectionFitTerm& setWeight(Real weight)
      {
        m_weight = std::max(Real(0), weight);
        return *this;
      }

      Real getWeight() const { return m_weight; }

      InterfaceNormalAdvectionFitTerm& setNormalization(Real cSigma)
      {
        m_normalization = (cSigma > Real(0)) ? cSigma : Real(1);
        return *this;
      }

      Real getNormalization() const { return m_normalization; }

      InterfaceNormalAdvectionFitTerm& setQuadratureOrder(size_t order)
      {
        m_quadratureOrder = order;
        return *this;
      }

      size_t getQuadratureOrder() const { return m_quadratureOrder; }

      template <class FES, class Data, class TestFES>
      auto residual(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TestFunction<TestFES>& v) const
      {
        return InterfaceNormalAdvectionFitResidualIntegrator<
            Variational::GridFunction<FES, Data>, TargetFunction>(
            u, v, m_target, m_interfaceAttribute,
            m_weight / m_normalization, m_quadratureOrder);
      }

      template <class FES, class Data, class Solution, class TrialFES, class TestFES>
      auto tangent(
          const Variational::GridFunction<FES, Data>& u,
          const Variational::TrialFunction<Solution, TrialFES>& du,
          const Variational::TestFunction<TestFES>& v) const
      {
        return InterfaceNormalAdvectionFitTangentIntegrator<
            Variational::GridFunction<FES, Data>, TargetFunction>(
            u, du, v, m_interfaceAttribute,
            m_weight / m_normalization, m_quadratureOrder);
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

        const auto& qf = QF::PolytopeQuadratureFormula::get(
            m_quadratureOrder, Geometry::Polytope::Type::Segment);

        Real value = 0;
        for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
        {
          if (!isActiveInterfaceEdge(mesh, e))
            continue;
          const auto edgeIterator = mesh.getPolytope(1, e);
          const auto& edge = *edgeIterator;
          const auto& fes  = u.getFiniteElementSpace();
          const auto& fe   = fes.getFiniteElement(edge.getDimension(), edge.getIndex());
          const size_t localSize = fe.getCount();
          const auto displacement = localDisplacementCoefficients(u, edge);

          for (size_t q = 0; q < qf.getSize(); ++q)
          {
            const auto& rc    = qf.getPoint(q);
            const auto   J    = referenceGeometryJacobian(edge, rc);
            const Real arcLen = std::hypot(J(0, 0), J(1, 0));
            if (arcLen <= Real(1e-14))
              continue;
            const Math::SpatialPoint n{ -J(1, 0) / arcLen, J(0, 0) / arcLen };

            Math::SpatialPoint x;
            edge.getTransformation().transform(x, rc);
            const Real tgt = m_target(x);

            Real uDotN = Real(0);
            for (size_t i = 0; i < localSize; ++i)
            {
              const auto phi = basisVectorValue(fe, i, rc);
              Real p = Real(0);
              for (size_t c = 0; c < sdim; ++c)
                p += phi[static_cast<std::uint8_t>(c)] * n[c];
              uDotN += displacement[i] * p;
            }
            const Real r = uDotN - tgt;
            value += qf.getWeight(q) * arcLen * r * r;
          }
        }
        return Real(0.5) * (m_weight / m_normalization) * value;
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

      TargetFunction m_target;
      Optional<Geometry::Attribute> m_interfaceAttribute;
      Real   m_weight          = 1;
      Real   m_normalization   = 1;
      size_t m_quadratureOrder = 4;
  };

  template <class TargetFunction>
  InterfaceNormalAdvectionFitTerm(
      TargetFunction,
      Optional<Geometry::Attribute>,
      Real) -> InterfaceNormalAdvectionFitTerm<TargetFunction>;

  template <class TargetFunction>
  InterfaceNormalAdvectionFitTerm(
      TargetFunction,
      Optional<Geometry::Attribute>) -> InterfaceNormalAdvectionFitTerm<TargetFunction>;

  template <class TargetFunction>
  InterfaceNormalAdvectionFitTerm(
      TargetFunction) -> InterfaceNormalAdvectionFitTerm<TargetFunction>;

}

#endif
