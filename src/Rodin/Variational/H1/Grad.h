/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_GRAD_H
#define RODIN_VARIATIONAL_H1_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for H1 (higher-order Lagrange) functions.
 *
 * For H1<K> functions, the gradient is polynomial of degree K-1 on each element:
 * @f[
 *   \nabla u|_K = \sum_{i=1}^{n} u_i \nabla \phi_i
 * @f]
 * where @f$ \phi_i @f$ are the H1<K> basis functions.
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  /**
   * @brief Base class for Grad classes.
   */
  template <class Operand, class Derived>
  class GradBase;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a GridFunction on H1<K> space
   */
  template <size_t K, class Scalar, class Mesh, class Data>
  class Grad<GridFunction<H1<K, Scalar, Mesh>, Data>> final
    : public GradBase<GridFunction<H1<K, Scalar, Mesh>, Data>, Grad<GridFunction<H1<K, Scalar, Mesh>, Data>>>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = GradBase<OperandType, Grad<OperandType>>;

      /**
       * @brief Constructs the gradient of an H1<K> function @f$ u @f$.
       * @param[in] u H1<K> GridFunction
       *
       * @note The gradient is a polynomial of degree K-1 on each element for H1<K> functions.
       */
      Grad(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Grad object to copy
       */
      Grad(const Grad& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Grad object to move from
       */
      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates the gradient at a given point.
       * @param[out] out Output spatial vector for gradient
       * @param[in] p Point at which to evaluate gradient
       *
       * Computes @f$ \nabla u(p) @f$ using H1<K> basis function gradients.
       * Handles evaluation on faces by projecting to adjacent cells.
       */
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
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
          else
          {
            assert(inc.size() == 2);
            const auto& traceDomain = this->getTraceDomain();
            assert(traceDomain.size() > 0);
            if (traceDomain.size() == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "No trace domain provided: "
                << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
                << ". Grad at an interface with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
                if (traceDomain.count(tracePolytope->getAttribute()))
                {
                  Math::SpatialPoint rc;
                  tracePolytope->getTransformation().inverse(rc, pc);
                  const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                  this->interpolate(out, np);
                  return;
                }
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          static thread_local SpatialVectorType s_res;

          s_res.resize(d);
          s_res.setZero();

          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto& basis = fe.getBasis(local);
            s_res += gf[fes.getGlobalIndex({d, i}, local)] * basis.getGradient()(rc);
          }
          out = p.getJacobianInverse().transpose() * s_res;
        }
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a ShapeFunction on H1<K> space
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, SpaceType>>>
  {
    public:
      /// Finite element space type
      using FESType = H1<K, Scalar, Mesh>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// Type of scalar values in the finite element space
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType = Math::Vector<ScalarType>;

      /// Operand type
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_p(nullptr)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_p(other.m_p),
          m_gradient(other.m_gradient)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_p(std::exchange(other.m_p, nullptr)),
          m_gradient(std::move(other.m_gradient))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      Grad& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;
        m_p = &p;
        const auto& polytope = p.getPolytope();
        const auto& rc = p.getReferenceCoordinates();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        m_gradient.resize(count);
        for (size_t local = 0; local < count; local++)
        {
          const auto& basis = fe.getBasis(local);
          m_gradient[local] = basis.getGradient()(rc);
        }
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return getPoint().getJacobianInverse().transpose() * m_gradient[local];
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      const Geometry::Point* m_p;
      std::vector<Math::SpatialVector<ScalarType>> m_gradient;
  };
}

#endif
