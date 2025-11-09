/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_POW_H
#define RODIN_VARIATIONAL_POW_H

/**
 * @file
 * @brief Power (exponentiation) operations for functions.
 */

#include <Rodin/Math/Common.h>

#include "ForwardDecls.h"

#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup PowSpecializations Pow Template Specializations
   * @brief Template specializations of the Pow class.
   * @see Pow
   */

  /**
   * @ingroup PowSpecializations
   * @brief Represents the power function.
   *
   * This class represents the exponentiation operation applied to a function:
   * @f[
   *  \text{Pow}(f, p)(x) = f(x)^p
   * @f]
   *
   * The base function @f$ f @f$ must be scalar-valued, and the exponent @f$ p @f$
   * is a compile-time constant. The result is always a real-valued scalar function.
   *
   * @tparam BaseDerived Type of the base function
   * @tparam Number Type of the exponent (must be arithmetic)
   *
   * @note Common use cases include quadratic functions (@f$ p = 2 @f$), cubic
   * functions (@f$ p = 3 @f$), and fractional powers for regularization.
   *
   * @see RealFunctionBase, Sqrt, Exp
   */
  template <class BaseDerived, class Number>
  class Pow<FunctionBase<BaseDerived>, Number> final
    : public RealFunctionBase<Pow<FunctionBase<BaseDerived>, Number>>
  {
    static_assert(std::is_arithmetic_v<Number>);
    public:
      /// Type of base
      using BaseType = FunctionBase<BaseDerived>;

      /// Type of exponent
      using ExponentType = Number;

      /// Parent class
      using Parent = RealFunctionBase<Pow<FunctionBase<BaseDerived>, Number>>;

      /**
       * @brief Constructs the power object
       * @param[in] s Base value
       * @param[in] p Power
       */
      constexpr
      Pow(const BaseType& s, ExponentType p)
        : m_s(s.copy()), m_p(p)
      {}

      constexpr
      Pow(const Pow& other)
        : Parent(other),
          m_s(other.m_s->copy()), m_p(other.m_p)
      {}

      constexpr
      Pow(Pow&& other)
        : Parent(std::move(other)),
          m_s(std::move(other.m_s)),
          m_p(std::move(other.m_p))
      {}

      template <class ... Args>
      constexpr
      Pow& traceOf(const Args& ... args)
      {
        m_s->traceOf(args...);
        return *this;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::pow(this->getBase().getValue(p), getExponent());
      }

      const BaseType& getBase() const
      {
        return *m_s;
      }

      const ExponentType& getExponent() const
      {
        return m_p;
      }

      Pow* copy() const noexcept override
      {
        return new Pow(*this);
      }

    private:
      std::unique_ptr<BaseType> m_s;
      const ExponentType m_p;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class BaseDerived, class Number>
  Pow(const FunctionBase<BaseDerived>&, Number) -> Pow<FunctionBase<BaseDerived>, Number>;

  template <class NestedDerived, class Number>
  auto pow(const FunctionBase<NestedDerived>& f, Number exponent)
  {
    return Pow(f, exponent);
  }
}

#endif
