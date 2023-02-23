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
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<Ts...>, Space>>, H1<Ts...>, Space>
  {
    public:
      using FES = H1<Ts...>;
      using Operand = ShapeFunction<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<Ts...>, Space>>, FES, Space>;

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
      size_t getDOFs(const Geometry::Simplex& simplex) const
      {
        return m_u.get().getDOFs(simplex);
      }

      auto getOperator(const Geometry::Point& p) const
      {
        return void();
        // Math::Vector div = fe.getGradient(p.getVector(Geometry::Point::Coordinates::Reference)).reshaped();
        // return TensorBasis(div);
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
