/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SINH_H
#define RODIN_VARIATIONAL_SINH_H

/**
 * @file Sinh.h
 * @brief Hyperbolic sine function operator for scalar functions.
 */

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup SinhSpecializations Sinh Template Specializations
   * @brief Template specializations of the Sinh class.
   * @see Sinh
   */

  /**
   * @brief Hyperbolic sine function operator for real-valued scalar functions.
   *
   * Applies the hyperbolic sine function pointwise to a given function:
   * @f[
   *    \text{Sinh}(f)(x) = \sinh(f(x)) = \frac{e^{f(x)} - e^{-f(x)}}{2}
   * @f]
   *
   * Common applications include:
   * - Solutions to hyperbolic PDEs
   * - Catenary curves and hanging cables
   * - Heat transfer in semi-infinite domains
   * - Exact solutions with exponential growth/decay
   *
   * @note Always returns a real-valued function for real input.
   * @see Cosh, Sin
   */

  /**
   * @ingroup SinhSpecializations
   */
  template <class NestedDerived>
  class Sinh<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Sinh<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Sinh<FunctionBase<NestedDerived>>>;

      /**
       * @brief Constructs hyperbolic sine operator from a function.
       * @param v Function to apply sinh to
       */
      Sinh(const OperandType& v)
        : m_operand(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Sinh operator to copy
       */
      Sinh(const Sinh& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Sinh operator to move
       */
      Sinh(Sinh&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Restricts sinh operator to a trace.
       * @param attr Attribute to restrict to
       * @return Reference to this sinh operator
       */
      constexpr
      Sinh& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      /**
       * @brief Restricts sinh operator to multiple trace attributes.
       * @param attrs Set of attributes to restrict to
       * @return Reference to this sinh operator
       */
      constexpr
      Sinh& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      /**
       * @brief Evaluates @f$ \sinh(f(p)) @f$ at a point.
       * @param p Point at which to evaluate
       * @return Hyperbolic sine of operand function value
       */
      Real getValue(const Geometry::Point& p) const
      {
        return Math::sinh(getOperand().getValue(p));
      }

      /**
       * @brief Gets the operand function.
       * @return Reference to operand function
       */
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Creates a polymorphic copy of this sinh operator.
       * @return Pointer to copy
       */
      Sinh* copy() const noexcept override
      {
        return new Sinh(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Sinh(const FunctionBase<NestedDerived>&) -> Sinh<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Sinh.
   */
  template <class NestedDerived>
  auto sinh(const FunctionBase<NestedDerived>& f)
  {
    return Sinh(f);
  }
}

#endif
