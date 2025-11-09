/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COSH_H
#define RODIN_VARIATIONAL_COSH_H

/**
 * @file Cosh.h
 * @brief Hyperbolic cosine function operator for scalar functions.
 */

#include "Rodin/Math/Common.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup CoshSpecializations Cosh Template Specializations
   * @brief Template specializations of the Cosh class.
   * @see Cosh
   */

  /**
   * @brief Hyperbolic cosine function operator for real-valued scalar functions.
   *
   * Applies the hyperbolic cosine function pointwise to a given function:
   * @f[
   *    \text{Cosh}(f)(x) = \cosh(f(x)) = \frac{e^{f(x)} + e^{-f(x)}}{2}
   * @f]
   *
   * Common applications include:
   * - Solutions to hyperbolic PDEs
   * - Catenary curves and hanging cables
   * - Temperature distributions in infinite domains
   * - Exact solutions with exponential behavior
   *
   * @note Always returns a real-valued function for real input.
   * @note Always positive: @f$ \cosh(x) \geq 1 @f$ for all real @f$ x @f$.
   * @see Sinh, Cos
   */

  /**
   * @ingroup CoshSpecializations
   */
  template <class NestedDerived>
  class Cosh<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Cosh<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Cosh<FunctionBase<NestedDerived>>>;

      Cosh(const OperandType& v)
        : m_operand(v.copy())
      {}

      Cosh(const Cosh& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Cosh(Cosh&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      template <class ... Args>
      Cosh& traceOf(const Args& ... args)
      {
        m_operand->traceOf(args...);
        return *this;
      }

      auto getValue(const Geometry::Point& p) const
      {
        return Math::cosh(m_operand->getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      Cosh* copy() const noexcept override
      {
        return new Cosh(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Cosh(const FunctionBase<NestedDerived>&) -> Cosh<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Cosh.
   */
  template <class NestedDerived>
  auto cosh(const FunctionBase<NestedDerived>& f)
  {
    return Cosh(f);
  }
}

#endif
