/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_GRADIENT_H
#define RODIN_VARIATIONAL_P1_GRADIENT_H

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Grad.h"

#include "GridFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 GridFunction
   */
  template <class ... Ts>
  class Grad<GridFunction<P1<Scalar, Ts...>>> final
    : public VectorFunctionBase<Grad<GridFunction<P1<Scalar, Ts...>>>>
  {
    public:
      /// Operand type
      using Operand = GridFunction<P1<Scalar, Ts...>>;

      /// Parent class
      using Parent = VectorFunctionBase<Grad<GridFunction<P1<Scalar, Ts...>>>>;

      /**
       * @brief Constructs the gradient of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      Grad(const Operand& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      inline
      constexpr
      Grad& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      Grad& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      auto getValue(const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& gf = m_u.get();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& rc = p.getReferenceCoordinates();
        Math::SpatialVector grad(d);
        Math::SpatialVector res(d);
        res.setZero();
        for (size_t local = 0; local < fe.getCount(); local++)
        {
          fe.getGradient(local)(grad, rc);
          res += gf.getValue(fes.getGlobalIndex({d, i}, local)) * grad;
        }
        return p.getJacobianInverse().transpose() * this->object(std::move(res));
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a P1 GridFunction
   */
  template <class ... Ts>
  Grad(const GridFunction<P1<Ts...>>&) -> Grad<GridFunction<P1<Ts...>>>;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 ShapeFunction
   */
  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>, P1<Scalar, Ps...>, Space>
  {
    public:
      /// Finite element space type
      using FES = P1<Scalar, Ps...>;

      /// Operand type
      using Operand = ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Grad<Operand>, FES, Space>;

      Grad(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getOperand().getFiniteElementSpace().getMesh().getSpaceDimension(), 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const auto& fes = this->getFiniteElementSpace();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t dofs = fe.getCount();
        const auto& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        m_gradient.resize(dofs);
        for (size_t local = 0; local < dofs; local++)
          m_gradient[local] = fe.getGradient(local)(rc);
        return TensorBasis(dofs,
            [&](size_t local) { return p.getJacobianInverse().transpose() * m_gradient[local]; });
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
      mutable std::vector<Math::SpatialVector> m_gradient;
  };

  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>;
}

#endif

