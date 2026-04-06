/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_DERIVATIVE_H
#define RODIN_VARIATIONAL_H1_DERIVATIVE_H

/**
 * @file Derivative.h
 * @brief Partial derivative operator specialization for H1 functions.
 *
 * This file provides specialized implementations of the partial derivative
 * operator for H1 GridFunctions. For an H1<K> function @f$ u @f$, the
 * partial derivative with respect to coordinate @f$ x_i @f$ is:
 * @f[
 *   \frac{\partial u}{\partial x_i} = (\nabla u) \cdot \mathbf{e}_i
 * @f]
 *
 * For H1<K> elements, partial derivatives are polynomial of degree K-1 on
 * each element.
 *
 * @see Derivative, Grad, H1
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Derivative.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "H1Element.h"

namespace Rodin::FormLanguage
{
  template <size_t K, class Scalar, class Mesh, class Data>
  struct Traits<Variational::Derivative<Variational::GridFunction<Variational::H1<K, Scalar, Mesh>, Data>>>
  {
    using FESType = Variational::H1<K, Scalar, Mesh>;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Scalar;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DerivativeSpecializations
   * @brief Derivative of an H1<K> GridFunction.
   *
   * Computes the partial derivative @f$ \frac{\partial u}{\partial x_i} @f$
   * of an H1<K> grid function @f$ u @f$. This is equivalent to the @f$ i @f$-th
   * component of the gradient @f$ \nabla u @f$.
   *
   * @tparam K Polynomial degree
   * @tparam Scalar Scalar value type
   * @tparam Mesh Mesh type
   * @tparam Data Data storage type
   */
  template <size_t K, class Scalar, class Mesh, class Data>
  class Derivative<GridFunction<H1<K, Scalar, Mesh>, Data>> final
    : public DerivativeBase<
        GridFunction<H1<K, Scalar, Mesh>, Data>,
        Derivative<GridFunction<H1<K, Scalar, Mesh>, Data>>>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;

      using RangeType = Scalar;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = DerivativeBase<OperandType, Derivative<OperandType>>;

      /**
       * @brief Constructs the partial derivative of an H1<K> function.
       * @param[in] i Coordinate index (0 = x, 1 = y, 2 = z)
       * @param[in] u H1 GridFunction
       */
      Derivative(size_t i, const OperandType& u)
        : Parent(u),
          m_i(i)
      {}

      /**
       * @brief Copy constructor.
       */
      Derivative(const Derivative& other)
        : Parent(other),
          m_i(other.m_i)
      {}

      /**
       * @brief Move constructor.
       */
      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_i(other.m_i)
      {}

      /**
       * @brief Interpolates the partial derivative at a given point.
       * @param[out] out Output scalar for derivative value
       * @param[in] p Point at which to evaluate
       *
       * Computes @f$ \frac{\partial u}{\partial x_i}(p) @f$ by computing
       * the full gradient and extracting the i-th component.
       * Handles evaluation on faces by projecting to adjacent cells.
       */
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();

        if (d == meshDim - 1) // Evaluating on a face
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc = conn.getIncidence({ meshDim - 1, meshDim }, i);
          const auto& pc = p.getPhysicalCoordinates();
          assert(inc.size() == 1 || inc.size() == 2);

          if (inc.size() == 1)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, *inc.begin());
            Math::SpatialPoint rc;
            tracePolytope->getTransformation().inverse(rc, pc);
            const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
            this->interpolate(out, np);
            return;
          }

          assert(inc.size() == 2);
          const auto& traceDomain = this->getTraceDomain();
          if (traceDomain.size() == 0)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "No trace domain provided: "
              << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
              << ". Derivative at an interface with no trace domain is undefined."
              << Alert::Raise;
          }

          for (auto& idx : inc)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
            const auto a = tracePolytope->getAttribute();
            if (a && traceDomain.count(*a))
            {
              Math::SpatialPoint rc;
              tracePolytope->getTransformation().inverse(rc, pc);
              const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
              this->interpolate(out, np);
              return;
            }
          }

          UndeterminedTraceDomainException(
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end())
            << Alert::Raise;
          return;
        }
        else // Evaluating on a cell
        {
          assert(d == mesh.getDimension());

          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();

          SpatialVectorType grad(d);
          grad.setZero();

          for (size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto& basis = fe.getBasis(local);
            grad += gf[fes.getGlobalIndex({d, i}, local)] * basis.getGradient()(rc);
          }

          SpatialVectorType transformed = p.getJacobianInverse().transpose() * grad;
          out = transformed(m_i);
        }
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = H1Element<K, ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      Derivative* copy() const noexcept override
      {
        return new Derivative(*this);
      }

    private:
      size_t m_i;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Derivative of an H1 GridFunction.
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  Derivative(size_t, const GridFunction<H1<K, Scalar, Mesh>, Data>&)
    -> Derivative<GridFunction<H1<K, Scalar, Mesh>, Data>>;
}

#endif
