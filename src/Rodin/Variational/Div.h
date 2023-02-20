#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"

#include "H1.h"
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
  class Div<ShapeFunction<NestedDerived, H1<Ts...>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<Ts...>, Space>>, Space>
  {
    public:
      using Operand = ShapeFunction<NestedDerived, H1<Ts...>, Space>;
      using Parent = ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<Ts...>, Space>>, Space>;

      /**
       * @brief Constructs Div object
       * @param[in] u ShapeFunction to be differentiated
       */
      constexpr
      Div(Operand& u)
        : m_u(u)
      {}

      constexpr
      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      constexpr
      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(other.m_u)
      {}

      inline
      constexpr
      const auto& getLeaf() const
      {
        return m_u.get().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { 1, 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return m_u.get().getDOFs(element);
      }

      auto getOperator(ShapeComputator& compute, const Geometry::Point& p) const
      {
        const auto& element = p.getSimplex();
        auto& trans = p.getSimplex().getTransformation();
        const auto& fe = getFiniteElementSpace().getFiniteElement(element);
        const auto& dshape = compute.getPhysicalDShape(fe, trans, trans.GetIntPoint());
        const size_t dofs = getDOFs(element);
        return TensorBasis(dofs, [](size_t i) -> Scalar { return dshape.GetData()[i]; });
      }

      inline
      constexpr
      auto& getFiniteElementSpace()
      {
        return m_u.get().getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return m_u.get().getFiniteElementSpace();
      }

    private:
      std::reference_wrapper<Operand> m_u;
  };

  template <class NestedDerived, ShapeFunctionSpaceType Space, class ... Ts>
  Div(ShapeFunction<NestedDerived, H1<Ts...>, Space>&)
    -> Div<ShapeFunction<NestedDerived, H1<Ts...>, Space>>;
}

#endif
