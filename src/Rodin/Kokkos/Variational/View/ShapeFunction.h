/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_VIEW_SHAPEFUNCTION_H
#define RODIN_KOKKOS_VARIATIONAL_VIEW_SHAPEFUNCTION_H

/**
 * @file
 * @brief Device evaluation and composition of trial and test functions.
 */

#include "Common.h"
#include "FiniteElementSpace.h"
#include "Function.h"
#include "GridFunction.h"

namespace Rodin::Kokkos::Variational
{
  /// @brief Forwards a shape function to its derived type.
  template <class Derived, class FES, ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::ShapeFunctionBase<Derived, FES, Space>>
    : View<Derived>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::ShapeFunctionBase<Derived, FES, Space>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Scalar basis of a device space.
  template <class Derived, ScalarSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::ShapeFunction<Derived, FES, Space>& node)
      : dofs(node.getFiniteElementSpace())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION Real getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component == 0);
      return dofs.getBasis(ip, dof);
    }

    /// @brief Device degrees of freedom.
    ::Rodin::Kokkos::View<FES> dofs;
  };

  /// @brief Physical gradient of the scalar basis of a device space.
  template <class Derived, ScalarSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Grad<
      ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>& node)
      : dofs(node.getOperand().getFiniteElementSpace())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Gradient.
     */
    KOKKOS_INLINE_FUNCTION decltype(auto) getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component == 0);
      return dofs.getGradient(ip, dof);
    }

    /// @brief Device degrees of freedom.
    ::Rodin::Kokkos::View<FES> dofs;
  };

  /// @brief Function multiplying a shape function on the left.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
    ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(
      const ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
        ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Scaled basis value.
     */
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return mult(evaluate(lhs, ip), rhs.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using ValueType = decltype(evaluate(lhs, ip));
      using BasisType = decltype(bindBasis(rhs, ip));
      return Internal::LeftScaledBasis<ValueType, BasisType>{
        evaluate(lhs, ip), bindBasis(rhs, ip)};
    }

    /// @brief Function factor.
    View<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Shape function factor.
    View<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /// @brief Function multiplying a shape function on the right.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::Mult<
    ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
    ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Mult<
      ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Scaled basis value.
     */
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return mult(lhs.getBasis(ip, dof, component), evaluate(rhs, ip));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using ValueType = decltype(evaluate(rhs, ip));
      using BasisType = decltype(bindBasis(lhs, ip));
      return Internal::RightScaledBasis<BasisType, ValueType>{
        bindBasis(lhs, ip), evaluate(rhs, ip)};
    }

    /// @brief Shape function factor.
    View<::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>> lhs;

    /// @brief Function factor.
    View<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Dot product of a function and a shape function.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
    ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Dot product value.
     */
    KOKKOS_INLINE_FUNCTION Real getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return dot(evaluate(lhs, ip), rhs.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using ValueType = decltype(evaluate(lhs, ip));
      using BasisType = decltype(bindBasis(rhs, ip));
      return Internal::DotBasis<ValueType, BasisType>{
        evaluate(lhs, ip), bindBasis(rhs, ip)};
    }

    /// @brief Function factor.
    View<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Shape function factor.
    View<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /**
   * @brief Vector basis of a device space.
   *
   * Local index @f$ a\,d + r @f$ is the basis function
   * @f$ \phi_a e_r @f$, as on the host.
   */
  template <class Derived, VectorSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::ShapeFunction<Derived, FES, Space>& node)
      : dofs(node.getFiniteElementSpace())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return dofs.getBasis(ip, dof, component);
    }

    /// @brief Device degrees of freedom.
    ::Rodin::Kokkos::View<FES> dofs;
  };

  /**
   * @brief Physical Jacobian of the vector basis of a device space.
   *
   * Row @f$ r @f$ of the Jacobian of @f$ \phi_a e_r @f$ is
   * @f$ \nabla \phi_a @f$, the other rows vanish.
   */
  template <class Derived, VectorSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::Jacobian<
    ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Jacobian<
      ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>& node)
      : dofs(node.getFiniteElementSpace())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Jacobian.
     */
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return dofs.getJacobian(ip, dof, component);
    }

    /// @brief Device degrees of freedom.
    ::Rodin::Kokkos::View<FES> dofs;
  };

  /**
   * @brief Divergence of the vector basis of a device space.
   *
   * The divergence of @f$ \phi_a e_r @f$ is @f$ \partial_r \phi_a @f$.
   */
  template <class Derived, VectorSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Div<
      ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>& node)
      : dofs(node.getFiniteElementSpace())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Divergence.
     */
    KOKKOS_INLINE_FUNCTION Real getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return dofs.getDivergence(ip, dof, component);
    }

    /// @brief Device degrees of freedom.
    ::Rodin::Kokkos::View<FES> dofs;
  };

  /// @brief Sum of two shape functions.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::Sum<
    ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
    ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Sum<
      ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return add(lhs.getBasis(ip, dof, component), rhs.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using LHSType = decltype(bindBasis(lhs, ip));
      using RHSType = decltype(bindBasis(rhs, ip));
      return Internal::SumBasis<LHSType, RHSType>{
        bindBasis(lhs, ip), bindBasis(rhs, ip)};
    }

    /// @brief Left operand.
    View<::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>> lhs;

    /// @brief Right operand.
    View<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /// @brief Negated shape function.
  template <class NestedDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::UnaryMinus<
    ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::UnaryMinus<
      ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>& node)
      : operand(node.getOperand())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return mult(Real(-1), operand.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using BasisType = decltype(bindBasis(operand, ip));
      return Internal::NegatedBasis<BasisType>{bindBasis(operand, ip)};
    }

    /// @brief Operand.
    View<::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>> operand;
  };

  /// @brief Transposed matrix valued shape function.
  template <class NestedDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct View<::Rodin::Variational::Transpose<
    ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    View(const ::Rodin::Variational::Transpose<
      ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>& node)
      : operand(node.getOperand())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION Geometry::Matrix getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return transpose(operand.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using BasisType = decltype(bindBasis(operand, ip));
      return Internal::TransposedBasis<BasisType>{bindBasis(operand, ip)};
    }

    /// @brief Operand.
    View<::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>> operand;
  };
}

#endif
