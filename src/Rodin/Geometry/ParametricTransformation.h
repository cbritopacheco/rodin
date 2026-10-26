/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_PARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_PARAMETRICTRANSFORMATION_H

/**
 * @file
 * @brief Parametric transformation for finite-element mesh geometry.
 */

#include <type_traits>
#include <limits>

#include "Rodin/Alert/MemberFunctionException.h"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PointCloud.h"

#include "PolytopeTransformation.h"

#include "ForwardDecls.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Geometry
{
  /**
   * @brief Parametric finite-element transformation for polytopes.
   *
   * Given geometry-node coordinates @f$P_i@f$ and scalar finite-element basis
   * functions @f$\phi_i@f$, the physical coordinate map is
   * @f[
   *   x(r) = \sum_i P_i \phi_i(r).
   * @f]
   *
   * The same class covers affine P1 geometry and curved Pk geometry.
   *
   * # Quadrature evaluation
   *
   * A bulk Jacobian evaluation borrows the element's existing reference
   * tabulation for the duration of the call, then contracts it with this
   * transformation's control points. It retains no table reference or new
   * cache. Elements without tabulation use the pointwise interface. The
   * scalar accumulation order is unchanged in either case.
   */
  template <class FE>
  class ParametricTransformation final : public PolytopeTransformation
  {
    static_assert(std::is_same_v<typename FE::RangeType, Real>,
      "Type of finite element must be scalar valued.");

    friend class boost::serialization::access;

    public:
      /// @brief Parent class type.
      using Parent = PolytopeTransformation;
      using Parent::transform;
      using Parent::jacobian;
      using Parent::inverse;

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element (the control-point count must match the element's
       * node count).
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(Geometry::PointCloud&& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(std::move(fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element.
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(const Geometry::PointCloud& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element.
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(Geometry::PointCloud&& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element.
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(const PointCloud& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(std::move(fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ParametricTransformation(const ParametricTransformation& other)
        : Parent(other),
          m_pm(other.m_pm),
          m_fe(other.m_fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      ParametricTransformation(ParametricTransformation&& other)
        : Parent(std::move(other)),
          m_pm(std::move(other.m_pm)),
          m_fe(std::move(other.m_fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      size_t getOrder() const override
      {
        return m_fe.getOrder();
      }

      /**
       * @brief Returns the element's factor degree when available, otherwise
       * the conservative total-degree bound of the transformation interface.
       * @returns The element's factor degree when available, otherwise the conservative total-degree bound of the transformation interface.
       */
      size_t getFactorOrder() const override
      {
        if constexpr (requires { m_fe.getFactorOrder(); })
          return m_fe.getFactorOrder();
        else
          return Parent::getFactorOrder();
      }

      void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const override
      {
        const size_t pdim = getPhysicalDimension();
        assert(static_cast<size_t>(rc.size()) == getReferenceDimension());
        pc.resize(pdim);
        pc.setZero();
        const size_t count = m_fe.getCount();
        for (size_t local = 0; local < count; local++)
        {
          assert(pc.size() == m_pm[local].size());
          pc += m_pm[local] * m_fe.getBasis(local)(rc);
        }
      }

      void jacobian(
        Math::SpatialMatrix<Real>& pc, const Math::SpatialPoint& rc) const override
      {
        const size_t rdim = getReferenceDimension();
        assert(static_cast<size_t>(rc.size()) == rdim);
        const size_t pdim = getPhysicalDimension();
        pc.resize(pdim, rdim);
        pc.setZero();
        const size_t count = m_fe.getCount();
        for (size_t local = 0; local < count; local++)
        {
          const auto& basis = m_fe.getBasis(local);
          for (size_t i = 0; i < rdim; i++)
          {
            const auto derivative = basis.template getDerivative<1>(i);
            // The reference derivative is shared by all physical components.
            const Real value = derivative(rc);
            // Expose the spatial storage bound to keep this a small loop.
            for (size_t j = 0; j < pdim && j < Math::SpatialMatrix<Real>::MaxSize; j++)
              pc(j, i) += m_pm(j, local) * value;
          }
        }
      }

      /**
       * @brief Computes quadrature Jacobians using the element's reference table.
       *
       * For control coordinates @f$P_{ja}@f$, each matrix entry is
       * @f$J_{ji}(\hat x_q)=\sum_a P_{ja}\partial_i\phi_a(\hat x_q)@f$.
       * Only reference derivatives are shared across cells; mapped matrices
       * remain specific to the current transformation. Generic scalar elements
       * without a table retain the default pointwise evaluation.
       * Square three-dimensional maps instantiate the same contraction with
       * compile-time spatial bounds. This removes repeated small-loop control
       * without changing the local-basis accumulation order, reference samples,
       * padding initialization or ownership. Other dimensions use runtime bounds;
       * no affine-map assumption or derived-metric cache is introduced.
       * The square-three-dimensional instance accumulates in a local fixed-capacity
       * matrix and commits it once per sample. Runtime-bound instances write to
       * retained output directly, avoiding an extra copy for smaller matrices.
       *
       * @param[out] jacobians Jacobian matrices in quadrature-point order.
       * @param[in] qf Reference formula identifying the element tabulation.
       */
      void jacobian(std::vector<Math::SpatialMatrix<Real>>& jacobians,
        const QF::QuadratureFormula& qf) const override
      {
        if constexpr (requires { m_fe.getTabulation(qf); })
        {
          const auto& table = m_fe.getTabulation(qf);
          const size_t rdim = getReferenceDimension(), pdim = getPhysicalDimension();
          assert(table.getDimension() == rdim);
          assert(table.getCount() == m_fe.getCount());
          jacobians.resize(table.getSize());
          const auto contract = [&](auto referenceDimension, auto physicalDimension) {
            constexpr bool localAccumulation =
              std::is_same_v<decltype(referenceDimension),
                std::integral_constant<size_t, 3>>;
            for (size_t qp = 0; qp < table.getSize(); ++qp)
            {
              decltype(auto) matrix = [&]() -> decltype(auto) {
                if constexpr (localAccumulation)
                  return Math::SpatialMatrix<Real>(physicalDimension, referenceDimension);
                else
                {
                  auto& output = jacobians[qp];
                  output.resize(physicalDimension, referenceDimension);
                  output.setZero();
                  return (output);
                }
              }();
              for (size_t local = 0; local < table.getCount(); ++local)
              {
                const auto gradient = table.getGradient(qp, local);
                for (size_t i = 0; i < referenceDimension; ++i)
                {
                  const Real value = gradient[i];
                  for (size_t j = 0;
                       j < physicalDimension && j < Math::SpatialMatrix<Real>::MaxSize;
                       ++j)
                  {
                    matrix(j, i) += m_pm(j, local) * value;
                  }
                }
              }
              if constexpr (localAccumulation)
                jacobians[qp] = matrix;
            }
          };
          if (rdim == 3 && pdim == 3)
          {
            contract(
              std::integral_constant<size_t, 3>{}, std::integral_constant<size_t, 3>{});
          }
          else
            contract(rdim, pdim);
        }
        else
          Parent::jacobian(jacobians, qf);
      }

      /**
       * @brief Inverts the finite-element geometry map at a physical point.
       * @param[out] rc Reference coordinates, assigned after convergence.
       * @param pc Physical point on the image of the reference element.
       *
       * # Architecture
       *
       * Affine simplices retain the direct inverse. Other parametric maps
       * require solving @f$F(r)=p@f$; one affine correction is insufficient
       * even for bilinear quadrilaterals. The map is assumed regular and
       * locally injective. This operation inverts a supplied point; reference
       * containment and mesh-wide element selection belong to the locator.
       *
       * # Implementation
       *
       * Newton corrections solve @f$J\delta=p-F(r)@f$. Embedded maps use
       * @f$(J^T J)\delta=J^T(p-F(r))@f$. A common Jacobian scale prevents
       * determinant overflow or underflow without changing the correction.
       * Backtracking requires a smaller physical residual. Both physical
       * residual and reference correction must reach roundoff-sized tolerances;
       * a tiny element alone therefore cannot authorize an inaccurate inverse.
       * Singular or nonconvergent maps raise a contextual exception.
       */
      void inverse(Math::SpatialPoint& rc, const Math::SpatialPoint& pc) const override
      {
        const auto norm = [](const Math::SpatialPoint& value) {
          Real result = 0;
          for (size_t i = 0; i < value.size(); ++i)
            result = std::hypot(result, value[i]);
          return result;
        };
        if (pc.size() != getPhysicalDimension())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Physical coordinates have an invalid extent." << Alert::Raise;
        }
        const Polytope::Traits traits(m_fe.getGeometry());
        if (getReferenceDimension() > 0 && traits.isSimplex() && getFactorOrder() <= 1)
        {
          Parent::inverse(rc, pc);
          return;
        }

        if (!std::isfinite(norm(pc)))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Physical coordinates must be finite." << Alert::Raise;
        }

        // Bounded local solve budgets; convergence is checked independently.
        constexpr size_t MaximumIterations = 32;
        constexpr size_t MaximumBacktracks = 16;
        constexpr Real RoundoffFactor = 64;
        const Real accuracy = RoundoffFactor * std::numeric_limits<Real>::epsilon();
        const Real physicalTolerance = accuracy * std::max(Real(1), norm(pc));
        Math::SpatialPoint reference = traits.getCentroid(), mapped;
        transform(mapped, reference);
        Math::SpatialPoint residual = pc - mapped;
        Real residualNorm = norm(residual);
        for (size_t iteration = 0; iteration < MaximumIterations; ++iteration)
        {
          if (residualNorm == 0)
          {
            rc = reference;
            return;
          }
          if (getReferenceDimension() == 0)
            break;
          Math::SpatialMatrix<Real> matrix;
          jacobian(matrix, reference);
          Real scale = 0;
          for (size_t i = 0; i < matrix.rows(); ++i)
          {
            for (size_t j = 0; j < matrix.cols(); ++j)
              scale = std::max(scale, std::abs(matrix(i, j)));
          }
          if (!(scale > 0) || !std::isfinite(scale))
            break;
          Math::SpatialPoint rhs(residual.size());
          for (size_t i = 0; i < matrix.rows(); ++i)
          {
            rhs[i] = residual[i] / scale;
            for (size_t j = 0; j < matrix.cols(); ++j)
              matrix(i, j) /= scale;
          }
          if (matrix.rows() != matrix.cols())
          {
            rhs = matrix.transpose() * rhs;
            matrix = matrix.transpose() * matrix;
          }
          const Real determinant = matrix.determinant();
          if (determinant == 0 || !std::isfinite(determinant))
            break;
          const Math::SpatialPoint step = matrix.solve(rhs);
          const Real stepNorm = norm(step);
          if (!std::isfinite(stepNorm))
            break;
          const Real referenceTolerance = accuracy * std::max(Real(1), norm(reference));
          if (residualNorm <= physicalTolerance && stepNorm <= referenceTolerance)
          {
            rc = reference;
            return;
          }
          Real length = 1;
          bool accepted = false;
          for (size_t backtrack = 0; backtrack < MaximumBacktracks; ++backtrack)
          {
            const Math::SpatialPoint candidate = reference + length * step;
            transform(mapped, candidate);
            Math::SpatialPoint candidateResidual = pc - mapped;
            const Real candidateNorm = norm(candidateResidual);
            if (candidateNorm < residualNorm)
            {
              reference = candidate;
              residual = candidateResidual;
              residualNorm = candidateNorm;
              accepted = true;
              break;
            }
            length *= Real(0.5);
          }
          if (!accepted)
            break;
        }
        Alert::MemberFunctionException(*this, __func__)
          << "Parametric inverse did not converge on a regular image point."
          << Alert::Raise;
      }

      /**
       * @brief Gets the control-point coordinate matrix.
       * @returns The control-point coordinate matrix.
       */
      const PointCloud& getPointMatrix() const
      {
        return m_pm;
      }

      /**
       * @brief Serializes the transformation (for boost::serialization).
       * @param ar Serialization archive.
       * @param version Boost.Serialization class version; unused by this implementation.
       */
      template <class Archive>
      void serialize(Archive& ar, [[maybe_unused]] const unsigned int version)
      {
        ar& boost::serialization::base_object<PolytopeTransformation>(*this);
        ar & m_pm;
        ar & m_fe;
      }

      ParametricTransformation* copy() const noexcept override
      {
        return new ParametricTransformation(*this);
      }

    protected:
      /**
       * @brief Maps a logical reference sample without retaining reference storage.
       *
       * When its existing table is present, @f$x_K(\hat x_q)=
       * \sum_a X_{K,a}\phi_a(\hat x_q)@f$ is accumulated in pointwise order.
       * Table eviction, another thread, and generic elements use direct
       * evaluation at the point's owned coordinates instead.
       *
       * @param[out] pc Physical coordinates.
       * @param rc Owned reference sample coordinates.
       * @param formula Persistent immutable reference formula.
       * @param qp Logical sample index established by point construction.
       */
      void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc,
        const QF::QuadratureFormula& formula, size_t qp) const override
      {
        if constexpr (requires { m_fe.findTabulation(formula); })
        {
          if (const auto found = m_fe.findTabulation(formula))
          {
            const auto& table = found->get();
            assert(table.getCount() == m_fe.getCount());
            assert(qp < table.getSize());
            pc.resize(getPhysicalDimension());
            pc.setZero();
            for (size_t local = 0; local < table.getCount(); ++local)
              pc += m_pm[local] * table.getBasis(qp, local);
            return;
          }
        }
        transform(pc, rc);
      }

    private:
      // Boost constructs this empty state before loading all serialized fields.
      ParametricTransformation()
        : Parent(0, 0)
      {}

      PointCloud m_pm;
      FE m_fe;
  };
}

#endif
