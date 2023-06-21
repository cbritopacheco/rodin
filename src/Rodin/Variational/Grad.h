/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "ForwardDecls.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  // /**
  //  * @ingroup GradSpecializations
  //  */
  // template <class ... Ts>
  // class Grad<GridFunction<H1<Scalar, Ts...>>> final
  //   : public VectorFunctionBase<Grad<GridFunction<H1<Scalar, Ts...>>>>
  // {
  //   public:
  //     using Operand = GridFunction<H1<Scalar, Ts...>>;
  //     using Parent = VectorFunctionBase<Grad<GridFunction<H1<Scalar, Ts...>>>>;

  //     /**
  //      * @brief Constructs the gradient of an @f$ H^1 @f$ function
  //      * @f$ u @f$.
  //      * @param[in] u Grid function to be differentiated
  //      */
  //     Grad(const Operand& u)
  //       : m_u(u)
  //     {
  //       assert(u.getFiniteElementSpace().getVectorDimension() == 1);
  //     }

  //     Grad(const Grad& other)
  //       : Parent(other),
  //         m_u(other.m_u)
  //     {}

  //     Grad(Grad&& other)
  //       : Parent(std::move(other)),
  //         m_u(std::move(other.m_u))
  //     {}

  //     inline
  //     constexpr
  //     size_t getDimension() const
  //     {
  //       return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
  //     }

  //     inline
  //     constexpr
  //     Grad& traceOf(Geometry::Attribute attr)
  //     {
  //       Parent::traceOf(attr);
  //       return *this;
  //     }

  //     inline
  //     constexpr
  //     Grad& traceOf(const std::set<Geometry::Attribute>& attrs)
  //     {
  //       Parent::traceOf(attrs);
  //       return *this;
  //     }

  //     Math::Vector getValue(const Geometry::Point& p) const
  //     {
  //     }

  //     inline
  //     constexpr
  //     const Operand& getOperand() const
  //     {
  //       return m_u.get();
  //     }

  //     inline Grad* copy() const noexcept override
  //     {
  //       return new Grad(*this);
  //     }

  //   private:
  //     std::reference_wrapper<const Operand> m_u;
  // };

  // template <class ... Ts>
  // Grad(const GridFunction<H1<Ts...>>&) -> Grad<GridFunction<H1<Ts...>>>;

  /**
   * @ingroup GradSpecializations
   */
  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>, P1<Scalar, Ps...>, Space>
  {
    public:
      using FES = P1<Scalar, Ps...>;
      using Operand = ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>;
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
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return TensorBasis(fe.getCount(),
            [&](size_t local)
            { return p.getJacobianInverse().transpose() * this->object(fe.getGradient(local)(rc)); });
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>;
}

#endif
