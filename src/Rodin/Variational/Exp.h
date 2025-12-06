/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EXP_H
#define RODIN_VARIATIONAL_EXP_H

/**
 * @file
 * @brief Exponential function operations.
 */

#include "Rodin/Math/Common.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ExpSpecializations Exp Template Specializations
   * @brief Template specializations of the Exp class.
   * @see Exp
   */

  /**
   * @ingroup ExpSpecializations
   * @brief Exponential function operation.
   *
   * This class represents the exponential operation applied to a scalar function:
   * @f[
   *    \text{Exp}(f)(x) = e^{f(x)}
   * @f]
   *
   * The input function must be scalar-valued (real or complex). The result is
   * always a real-valued scalar function when the input is real.
   *
   * Common applications include:
   * - Growth and decay models
   * - Softmax functions in optimization
   * - Barrier functions in interior point methods
   *
   * @tparam NestedDerived Type of the operand function
   *
   * @see RealFunctionBase, Pow
   */
  template <class NestedDerived>
  class Exp<FunctionBase<NestedDerived>>
    : public RealFunctionBase<Exp<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Exp<OperandType>>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      /**
       * @brief Constructs exponential of a function.
       * @param v Function to exponentiate
       */
      Exp(const OperandType& v)
        : m_v(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Exp object to copy from
       */
      Exp(const Exp& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Exp object to move from
       */
      Exp(Exp&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      /**
       * @brief Restricts operand to a trace domain.
       * @param args Arguments for trace restriction
       * @returns Reference to this object
       */
      template <class ... Args>
      constexpr
      Exp& traceOf(const Args& ... args)
      {
        m_v->traceOf(args...);
        return *this;
      }

      /**
       * @brief Evaluates exponential at a point.
       * @param p Point at which to evaluate
       * @returns @f$ e^{f(p)} @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::exp(this->getOperand().getValue(p));
      }

      /**
       * @brief Gets the operand function.
       * @returns Reference to the function in the exponent
       */
      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Exp* copy() const noexcept override
      {
        return new Exp(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Exp(const FunctionBase<NestedDerived>&) -> Exp<FunctionBase<NestedDerived>>;

  template <class NestedDerived>
  constexpr auto
  exp(const FunctionBase<NestedDerived>& op)
  {
    return Exp(op);
  }
}

#endif



