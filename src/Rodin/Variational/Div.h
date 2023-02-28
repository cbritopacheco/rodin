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
  class Div<ShapeFunction<NestedDerived, H1<Math::Vector, Ts...>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<Math::Vector, Ts...>, Space>>, H1<Math::Vector, Ts...>, Space>
  {
    public:
      using FES = H1<Math::Vector, Ts...>;
      using Operand = ShapeFunction<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Div<ShapeFunction<NestedDerived, FES, Space>>, FES, Space>;

      /**
       * @brief Constructs Div object
       * @param[in] u ShapeFunction to be differentiated
       */
      constexpr
      Div(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      constexpr
      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      constexpr
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
      size_t getDOFs(const Geometry::Simplex& simplex) const
      {
        return getOperand().getDOFs(simplex);
      }

      inline
      TensorBasis<Scalar> getTensorBasis(const Geometry::Point& p) const
      {
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(p.getSimplex());
        const auto& inv = p.getJacobianInverse();
        const Eigen::TensorMap<const Eigen::Tensor<Scalar, 2>> lift(inv.data(), inv.rows(), inv.cols());
        static constexpr const Eigen::array<Eigen::IndexPair<int>, 1> dims = { Eigen::IndexPair<int>(2, 0) };
        Math::Tensor<1> basis = fe.getJacobian(p.getReference()).contract(lift, dims).trace(Eigen::array<int, 2>{1, 2});
        return Eigen::Map<Math::Vector>(basis.data(), basis.size());
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
  Div(const ShapeFunction<NestedDerived, H1<Math::Vector, Ts...>, Space>&)
    -> Div<ShapeFunction<NestedDerived, H1<Math::Vector, Ts...>, Space>>;
}

#endif
