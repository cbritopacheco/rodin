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
       * @brief Constructs the power operation.
       * @param s Base function @f$ f @f$
       * @param p Exponent @f$ p @f$
       */
      constexpr
      Pow(const BaseType& s, ExponentType p)
        : m_s(s.copy()), m_p(p)
      {}

      /**
       * @brief Copy constructor.
       * @param other Pow object to copy from
       */
      constexpr
      Pow(const Pow& other)
        : Parent(other),
          m_s(other.m_s->copy()), m_p(other.m_p)
      {}

      /**
       * @brief Move constructor.
       * @param other Pow object to move from
       */
      constexpr
      Pow(Pow&& other)
        : Parent(std::move(other)),
          m_s(std::move(other.m_s)),
          m_p(std::move(other.m_p))
      {}

      /**
       * @brief Restricts base to a trace domain.
       * @param args Arguments for trace restriction
       * @returns Reference to this object
       */
      template <class ... Args>
      constexpr
      Pow& traceOf(const Args& ... args)
      {
        m_s->traceOf(args...);
        return *this;
      }

      /**
       * @brief Evaluates power at a point.
       * @param p Point at which to evaluate
       * @returns @f$ f(p)^{\text{exponent}} @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::pow(this->getBase().getValue(p), getExponent());
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const auto o = getBase().getOrder(polytope);
        if (!o)
          return std::nullopt;
        if constexpr (std::is_integral_v<ExponentType>)
        {
          const auto rawExp = getExponent();
          if constexpr (std::is_signed_v<ExponentType>)
          {
            if (rawExp < 0)
              return std::nullopt;
          }
          if (rawExp == 0)
            return size_t{0};
          const auto exp = static_cast<size_t>(rawExp);
          // At this point exp >= 1 because the rawExp == 0 case has returned above.
          const auto limit = std::numeric_limits<size_t>::max() / exp;
          if (*o > limit)
            return std::nullopt;
          return (*o) * exp;
        }
        return std::nullopt;
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
