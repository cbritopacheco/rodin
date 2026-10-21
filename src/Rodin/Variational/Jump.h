/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jump.h
 * @brief Jump operator for interface discontinuities in DG methods.
 *
 * This file defines the Jump class, which computes the jump (difference) in
 * function values across element interfaces. This operator is fundamental for
 * Discontinuous Galerkin (DG) formulations where functions are allowed to be
 * discontinuous at element boundaries.
 */
#ifndef RODIN_VARIATIONAL_JUMP_H
#define RODIN_VARIATIONAL_JUMP_H

#include <algorithm>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"

#include "ForwardDecls.h"
#include "IntegrationPoint.h"
#include "ShapeFunction.h"
#include "InterfaceShapeFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup JumpSpecializations Jump Template Specializations
   * @brief Template specializations of the Jump class.
   * @see @ref Jump
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Jump "Jump<FunctionBase<FunctionDerived>>" | Jump of a traceable function across an interior interface. |
   * | @ref Jump "Jump<ShapeFunctionBase<NestedDerived, FES, Space>>" | Jump of a test or trial shape function across an interior interface. |
   */

  /**
   * @ingroup JumpSpecializations
   * @brief Jump of a FunctionBase instance across an interface.
   *
   * The Jump operator computes the jump (discontinuity) in a function's value
   * across an interface (face) between two adjacent elements:
   * @f[
   *   [\![u]\!] = u^+ - u^-
   * @f]
   * where @f$ u^+ @f$ and @f$ u^- @f$ denote the values from the two sides
   * of the interface (with a chosen orientation).
   *
   * ## Mathematical Foundation
   * In Discontinuous Galerkin (DG) methods, functions are allowed to be
   * discontinuous across element interfaces. The jump operator captures this
   * discontinuity and is used in:
   * - **Interior Penalty terms**: @f$ \int_{\Gamma_h} [\![u]\!] \cdot [\![v]\!] \, ds @f$
   * - **Consistency terms**: @f$ \int_{\Gamma_h} \{\!\{\nabla u\}\!\} \cdot [\![v]\!] \, ds @f$
   * - **Flux terms**: @f$ \int_{\Gamma_h} [\![u]\!] \cdot \{\!\{\nabla v\}\!\} \, ds @f$
   *
   * where @f$ \Gamma_h @f$ denotes the set of interior faces.
   *
   * ## Usage Example
   * ```cpp
   * // Evaluate jump of a function on an interface
   * auto jump_u = Jump(u);
   * ```
   *
   * @tparam FunctionDerived Type of the function being jumped
   *
   * @see @ref Average
   */
  template <class FunctionDerived>
  class Jump<FunctionBase<FunctionDerived>> final
    : public FunctionBase<Jump<FunctionBase<FunctionDerived>>>
  {
    public:
      /// @brief Type of the operand function
      using OperandType = FunctionBase<FunctionDerived>;

      /// @brief Parent class type
      using Parent = FunctionBase<Jump<FunctionBase<FunctionDerived>>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs a jump operator for the given function.
       * @param[in] op Function to compute the jump of across interfaces
       */
      constexpr Jump(const OperandType& op)
        : m_operand(op.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Jump operator to copy
       */
      constexpr Jump(const Jump& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Jump operator to move
       */
      constexpr Jump(Jump&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Gets the operand function.
       * @returns Const reference to the operand
       */
      constexpr const auto& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Evaluates the jump at a point on an interface.
       * @param[in] p Point on an interior face
       * @returns Jump value @f$ u^+(p) - u^-(p) @f$
       *
       * The point must be on an interior face with exactly two adjacent elements.
       */
      auto getValue(const Geometry::Point& p) const
      {
        assert(p.getPolytope().isFace());
        const auto& face = p.getPolytope();
        const size_t d = face.getDimension();
        const auto& mesh = face.getMesh();
        const auto& inc =
          mesh.getConnectivity().getIncidence({d, d + 1}, face.getIndex());
        assert(inc.size() == 2);
        const Index idx1 = *inc.begin();
        const Index idx2 = *std::next(inc.begin());
        const auto it1 = mesh.getPolytope(d + 1, idx1);
        const auto it2 = mesh.getPolytope(d + 1, idx2);
        const auto& pc = p.getPhysicalCoordinates();
        Math::SpatialPoint rc1;
        Math::SpatialPoint rc2;
        it1->getTransformation().inverse(rc1, pc);
        it2->getTransformation().inverse(rc2, pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
        const auto v1 = getOperand().getValue(p1);
        const auto v2 = getOperand().getValue(p2);
        if constexpr (requires { (v1 - v2).eval(); })
          return (v1 - v2).eval();
        else
          return v1 - v2;
      }

      /**
       * @brief Evaluates the expression at an integration point.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      auto getValue(const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        assert(p.getPolytope().isFace());
        const auto& face = p.getPolytope();
        const size_t d = face.getDimension();
        const auto& mesh = face.getMesh();
        const auto& inc =
          mesh.getConnectivity().getIncidence({d, d + 1}, face.getIndex());
        assert(inc.size() == 2);
        const Index idx1 = *inc.begin();
        const Index idx2 = *std::next(inc.begin());
        const auto it1 = mesh.getPolytope(d + 1, idx1);
        const auto it2 = mesh.getPolytope(d + 1, idx2);
        const auto& pc = p.getPhysicalCoordinates();
        Math::SpatialPoint rc1;
        Math::SpatialPoint rc2;
        it1->getTransformation().inverse(rc1, pc);
        it2->getTransformation().inverse(rc2, pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
        const auto* qf = ip.getQuadratureFormula();
        assert(qf);
        const IntegrationPoint ip1(p1, qf, ip.getIndex());
        const IntegrationPoint ip2(p2, qf, ip.getIndex());
        const auto v1 = getOperand().getValue(ip1);
        const auto v2 = getOperand().getValue(ip2);
        if constexpr (requires { (v1 - v2).eval(); })
          return (v1 - v2).eval();
        else
          return v1 - v2;
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param p Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(const Geometry::Polytope& p) const noexcept
      {
        return this->getOperand().getOrder(p);
      }

      /**
       * @brief Creates a copy of this jump operator.
       * @returns Pointer to newly allocated copy
       */
      Jump* copy() const noexcept override
      {
        return new Jump(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for Jump with FunctionBase.
   * @param op Function to compute the jump of across interfaces
   */
  template <class Derived>
  Jump(const FunctionBase<Derived>& op) -> Jump<FunctionBase<Derived>>;

  /**
   * @ingroup JumpSpecializations
   * @brief Jump of a ShapeFunctionBase instance across an interface.
   *
   * Computes the jump of shape (trial or test) functions across element
   * interfaces in finite element formulations. This is used in DG bilinear
   * and linear form assembly.
   *
   * ## Usage Example
   * ```cpp
   * // DG interior penalty term
   * auto penalty = InterfaceIntegral(Jump(u), Jump(v));
   * ```
   *
   * @tparam NestedDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Jump<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public InterfaceShapeFunctionBase<
        Jump<ShapeFunctionBase<NestedDerived, FES, Space>>,
        ShapeFunctionBase<NestedDerived, FES, Space>>
  {
    public:
      /// @brief Cell-local shape expression.
      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      /// @brief Two-sided evaluation base.
      using Parent = InterfaceShapeFunctionBase<Jump, OperandType>;

      /**
       * @brief Constructs a two-sided jump.
       * @param operand Cell-local shape expression.
       */
      explicit Jump(const OperandType& operand)
        : Parent(operand, 1, -1)
      {}

      /// @returns An independent expression copy, without borrowed evaluation state.
      Jump* copy() const noexcept override
      {
        return new Jump(*this);
      }
  };

  /**
   * @brief Deduction guide for Jump with ShapeFunctionBase.
   * @param op Function to compute the jump of across interfaces
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Jump(const ShapeFunctionBase<NestedDerived, FES, Space>& op)
    -> Jump<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
