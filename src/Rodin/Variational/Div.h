#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"
#include "Trace.h"

namespace Rodin::Variational
{
  /**
    * @defgroup DivSpecializations Div Template Specializations
    * @brief Template specializations of the Div class.
    * @see Div
    */

  /**
   * @ingroup DivSpecializations
   */
  template <class NestedDerived, ShapeFunctionSpaceType Space, class ... Ts>
  class Div<ShapeFunction<NestedDerived, P1<Math::Vector, Ts...>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, P1<Math::Vector, Ts...>, Space>>, P1<Math::Vector, Ts...>, Space>
  {
    public:
      using FES = P1<Math::Vector, Ts...>;
      using Operand = ShapeFunction<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Div<ShapeFunction<NestedDerived, FES, Space>>, FES, Space>;

      /**
       * @brief Constructs Div object
       * @param[in] u ShapeFunction to be differentiated
       */
      Div(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Div(Div&& other)
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
        return { 1, 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return getOperand().getDOFs(polytope);
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return TensorBasis<Scalar>(fe.getCount(),
            [&](size_t local) -> Scalar
            { return (p.getJacobianInverse().transpose() * fe.getJacobian(local)(rc)).trace(); });
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  template <class NestedDerived, ShapeFunctionSpaceType Space, class ... Ts>
  Div(const ShapeFunction<NestedDerived, P1<Math::Vector, Ts...>, Space>&)
    -> Div<ShapeFunction<NestedDerived, P1<Math::Vector, Ts...>, Space>>;
}

#endif
