/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TANGENT_H
#define RODIN_VARIATIONAL_TANGENT_H

/**
 * @file Tangent.h
 * @brief Tangent trigonometric function operator for scalar functions.
 */

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup TanSpecializations Tan Template Specializations
   * @brief Template specializations of the Tan class.
   * @see Tan
   */

  /**
   * @brief Tangent function operator for real-valued scalar functions.
   *
   * Applies the tangent function pointwise to a given function:
   * @f[
   *    \text{Tan}(f)(x) = \tan(f(x)) = \frac{\sin(f(x))}{\cos(f(x))}
   * @f]
   *
   * Common applications include:
   * - Trigonometric identities in manufactured solutions
   * - Nonlinear trigonometric problems
   * - Angle computations
   *
   * @note Undefined when @f$ \cos(f(x)) = 0 @f$ (at odd multiples of @f$ \pi/2 @f$).
   * @note Always returns a real-valued function for real input.
   * @see Sin, Cos
   */

  /**
   * @ingroup TanSpecializations
   */
  template <class NestedDerived>
  class Tan<FunctionBase<NestedDerived>>
    : public RealFunctionBase<Tan<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Tan<FunctionBase<NestedDerived>>>;

      /**
       * @brief Constructs tangent operator from a function.
       * @param v Function to apply tangent to
       */
      constexpr
      Tan(const OperandType& v)
        : m_operand(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Tangent operator to copy
       */
      constexpr
      Tan(const Tan& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Tangent operator to move
       */
      constexpr
      Tan(Tan&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Restricts tangent operator to a trace.
       * @param args Trace restriction arguments
       * @return Reference to this tangent operator
       */
      template <class ... Args>
      constexpr
      Tan& traceOf(Args&& ... args)
      {
        Parent::traceOf(std::forward<Args>(args)...);
        m_operand->traceOf(std::forward<Args>(args)...);
        return *this;
      }

      /**
       * @brief Evaluates @f$ \tan(f(p)) @f$ at a point.
       * @param p Point at which to evaluate
       * @return Tangent of operand function value
       */
      Real getValue(const Geometry::Point& p) const
      {
        return Math::tan(getOperand().getValue(p));
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
       * @brief Creates a polymorphic copy of this tangent operator.
       * @return Pointer to copy
       */
      Tan* copy() const noexcept
      override
      {
        return new Tan(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Tan(const FunctionBase<NestedDerived>&) -> Tan<FunctionBase<NestedDerived>>;

  template <class NestedDerived>
  auto tan(const FunctionBase<NestedDerived>& f)
  {
    return Tan(f);
  }
}

#endif

