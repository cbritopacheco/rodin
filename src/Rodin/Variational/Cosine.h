/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COS_H
#define RODIN_VARIATIONAL_COS_H

/**
 * @file Cosine.h
 * @brief Cosine trigonometric function operator for scalar functions.
 */

#include "Rodin/Math/Common.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup CosSpecializations Cos Template Specializations
   * @brief Template specializations of the Cos class.
   * @see Cos
   */

  /**
   * @brief Cosine function operator for real-valued scalar functions.
   *
   * Applies the cosine function pointwise to a given function:
   * @f[
   *    \text{Cos}(f)(x) = \cos(f(x))
   * @f]
   *
   * Common applications include:
   * - Periodic boundary conditions
   * - Wave equations and harmonic oscillators
   * - Fourier series representations
   * - Trigonometric manufactured solutions
   *
   * @note Always returns a real-valued function for real input.
   * @see Sin, Tan, Cosh
   */

  /**
   * @ingroup CosSpecializations
   */
  template <class NestedDerived>
  class Cos<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Cos<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Cos<FunctionBase<NestedDerived>>>;

      /**
       * @brief Constructs cosine operator for a function.
       * @param v Function to apply cosine to
       */
      Cos(const OperandType& v)
        : m_operand(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Cosine operator to copy
       */
      Cos(const Cos& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Cosine operator to move from
       */
      Cos(Cos&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Restricts evaluation to specified mesh attributes.
       * @param args Attribute arguments for trace restriction
       * @returns Reference to this object
       */
      template <class ... Args>
      constexpr
      Cos& traceOf(const Args& ... args)
      {
        m_operand->traceOf(args...);
        return *this;
      }

      /**
       * @brief Evaluates cosine at a point.
       * @param p Point at which to evaluate
       * @returns @f$ \cos(f(p)) @f$
       */
      Real getValue(const Geometry::Point& p) const
      {
        return Math::cos(getOperand().getValue(p));
      }

      /**
       * @brief Gets the operand function.
       * @returns Reference to the function @f$ f @f$
       */
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Polymorphic copy.
       * @returns Pointer to a copy of this object
       */
      Cos* copy() const noexcept override
      {
        return new Cos(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Cos(const FunctionBase<NestedDerived>&) -> Cos<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Cos.
   */
  template <class NestedDerived>
  auto cos(const FunctionBase<NestedDerived>& f)
  {
    return Cos(f);
  }
}

#endif
