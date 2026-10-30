/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Average.h
 * @brief Average operator for interface jump terms in DG methods.
 *
 * This file defines the Average class, which computes the average value of
 * a function across element interfaces. This operator is essential for
 * Discontinuous Galerkin (DG) formulations and interior penalty methods.
 */
#ifndef RODIN_VARIATIONAL_AVERAGE_H
#define RODIN_VARIATIONAL_AVERAGE_H

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
   * @defgroup AverageSpecializations Average Template Specializations
   * @brief Template specializations of the Average class.
   * @see @ref Average
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Average "Average<FunctionBase<FunctionDerived>>" | Arithmetic average of a traceable function across an interior interface. |
   * | @ref Average "Average<ShapeFunctionBase<NestedDerived, FES, Space>>" | Arithmetic average of a test or trial shape function across an interior interface. |
   */

  /**
   * @ingroup AverageSpecializations
   * @brief Average of a FunctionBase instance across an interface.
   *
   * The Average operator computes the average value of a function across an
   * interface (face) between two adjacent elements.
   * @f[
   *   \{\!\{u\}\!\} = \frac{1}{2}(u^+ + u^-)
   * @f]
   * where @f$ u^+ @f$ and @f$ u^- @f$ denote the values from the two sides
   * of the interface.
   *
   * # Mathematical Foundation
   * In DG methods, functions are allowed to be discontinuous across element
   * interfaces. The average operator is used to define consistent numerical
   * fluxes and penalization terms:
   * - **Interior Penalty**: @f$ \int_{\Gamma_h} \{\!\{\nabla u\}\!\} \cdot [\![v]\!] \, ds @f$
   * - **Consistency**: @f$ \int_{\Gamma_h} \{\!\{\nabla v\}\!\} \cdot [\![u]\!] \, ds @f$
   *
   * where @f$ \Gamma_h @f$ denotes the set of interior faces.
   *
   * # Usage Example
   * ```cpp
   * // DG interior penalty term
   * auto penalty = InterfaceIntegral(Average(Grad(u)), Jump(v));
   * ```
   *
   * @tparam FunctionDerived Type of the function being averaged
   *
   * @see @ref Jump
   */
  template <class FunctionDerived>
  class Average<FunctionBase<FunctionDerived>> final
    : public FunctionBase<Average<FunctionBase<FunctionDerived>>>
  {
    public:
      /// @brief Type of the operand function
      using OperandType = FunctionBase<FunctionDerived>;

      /// @brief Parent class type
      using Parent = FunctionBase<Average<FunctionBase<FunctionDerived>>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs an average operator for the given function.
       * @param[in] op Function to average across interfaces
       */
      constexpr Average(const OperandType& op)
        : m_operand(op.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Average operator to copy
       */
      constexpr Average(const Average& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Average operator to move
       */
      constexpr Average(Average&& other)
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
       * @brief Evaluates the average at a point on an interface.
       * @param[in] p Point on an interior face
       * @returns Average value @f$ \frac{1}{2}(u^+ + u^-) @f$
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
        if constexpr (requires { (0.5 * (v1 + v2)).eval(); })
          return (0.5 * (v1 + v2)).eval();
        else
          return 0.5 * (v1 + v2);
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
        const IntegrationPoint ip1(p1, qf, ip.getIndex());
        const IntegrationPoint ip2(p2, qf, ip.getIndex());
        const auto v1 = getOperand().getValue(ip1);
        const auto v2 = getOperand().getValue(ip2);
        if constexpr (requires { (0.5 * (v1 + v2)).eval(); })
          return (0.5 * (v1 + v2)).eval();
        else
          return 0.5 * (v1 + v2);
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
       * @brief Creates a copy of this average operator.
       * @returns Pointer to newly allocated copy
       */
      Average* copy() const noexcept override
      {
        return new Average(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for Average with FunctionBase.
   * @param op Shape function to average
   */
  template <class Derived>
  Average(const FunctionBase<Derived>& op) -> Average<FunctionBase<Derived>>;

  /**
   * @ingroup AverageSpecializations
   * @brief Average of a ShapeFunctionBase instance across an interface.
   *
   * Computes the average of shape (trial or test) functions across element
   * interfaces in finite element formulations. This is used in DG bilinear
   * and linear form assembly.
   *
   * # Usage Example
   * ```cpp
   * // DG consistency term
   * auto consistency = InterfaceIntegral(Average(Grad(u)), Jump(v));
   * ```
   *
   * @tparam NestedDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Average<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public InterfaceShapeFunctionBase<
        Average<ShapeFunctionBase<NestedDerived, FES, Space>>,
        ShapeFunctionBase<NestedDerived, FES, Space>>
  {
    public:
      /// @brief Cell-local shape expression.
      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      /// @brief Two-sided evaluation base.
      using Parent = InterfaceShapeFunctionBase<Average, OperandType>;

      /**
       * @brief Constructs a two-sided average.
       * @param operand Cell-local shape expression.
       */
      explicit Average(const OperandType& operand)
        : Parent(operand, 0.5, 0.5)
      {}

      /// @returns An independent expression copy, without borrowed evaluation state.
      Average* copy() const noexcept override
      {
        return new Average(*this);
      }
  };

  /**
   * @brief Deduction guide for Average with ShapeFunctionBase.
   * @param op Shape function to average
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Average(const ShapeFunctionBase<NestedDerived, FES, Space>& op)
    -> Average<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
