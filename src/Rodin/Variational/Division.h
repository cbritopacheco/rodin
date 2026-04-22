/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIVISION_H
#define RODIN_VARIATIONAL_DIVISION_H

/**
 * @file
 * @brief Division operation for functions.
 */

#include "ForwardDecls.h"
#include "Function.h"

namespace Rodin::Variational
{
  /**
   * @defgroup DivisionSpecializations Division Template Specializations
   * @brief Template specializations of the Division class.
   * @see Division
   */

  /**
   * @ingroup DivisionSpecializations
   * @brief Division of a function by another function.
   *
   * Represents the mathematical expression:
   * @f[
   *    \left(\frac{f}{g}\right)(x) = \frac{f(x)}{g(x)}
   * @f]
   * where @f$ f @f$ and @f$ g @f$ are functions with compatible ranges.
   *
   * @tparam LHSDerived Type of the numerator function
   * @tparam RHSDerived Type of the denominator function
   *
   * @note The denominator function @f$ g @f$ must not evaluate to zero on the
   * domain of interest.
   */
  template <class LHSDerived, class RHSDerived>
  class Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>
    : public FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent =
        FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      /**
       * @brief Constructs a division from numerator and denominator functions.
       * @param lhs Numerator function @f$ f @f$
       * @param rhs Denominator function @f$ g @f$
       */
      Division(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Division object to copy
       */
      Division(const Division& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Division object to move from
       */
      Division(Division&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Restricts both functions to the trace of a mesh object.
       * @param args Arguments specifying the trace
       * @returns Reference to this object
       */
      template <class ... Args>
      Division& traceOf(const Args& ... args)
      {
        m_lhs->traceOf(args...);
        m_rhs->traceOf(args...);
        return *this;
      }

      /**
       * @brief Gets the numerator function.
       * @returns Reference to numerator function @f$ f @f$
       */
      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      /**
       * @brief Gets the denominator function.
       * @returns Reference to denominator function @f$ g @f$
       */
      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      /**
       * @brief Evaluates the division at a point.
       * @param p Point at which to evaluate
       * @returns Value @f$ \frac{f(x)}{g(x)} @f$ at point @f$ x @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getLHS().getValue(p)) / this->object(getRHS().getValue(p));
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const auto lo = getLHS().getOrder(polytope);
        const auto ro = getRHS().getOrder(polytope);
        // Only polynomial if denominator is polynomial of order 0 (constant)
        if (!lo || !ro || *ro != 0)
          return std::nullopt;
        return lo;
      }

      /**
       * @brief Creates a copy of the division operation.
       * @returns Pointer to copied object
       */
      Division* copy() const noexcept final override
      {
        return new Division(*this);
      }

    private:
      std::unique_ptr<FunctionBase<LHSDerived>> m_lhs;
      std::unique_ptr<FunctionBase<RHSDerived>> m_rhs;
  };
  /**
   * @brief Deduction guide for Division.
   */
  template <class LHSDerived, class RHSDerived>
  Division(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @brief Division operator for two functions.
   * @param lhs Numerator function
   * @param rhs Denominator function
   * @returns Division object representing @f$ \frac{f}{g} @f$
   */
  template <class LHSDerived, class RHSDerived>
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(lhs, rhs);
  }

  /**
   * @brief Division of a function by a number.
   * @param lhs Function to divide
   * @param rhs Number to divide by
   * @returns Division object representing @f$ \frac{f}{c} @f$
   */
  template <class LHSDerived, class Number,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<LHSDerived, RealFunction<Number>>>>
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Division(lhs, RealFunction(rhs));
  }

  /**
   * @brief Division of a number by a function.
   * @param lhs Number in numerator
   * @param rhs Function in denominator
   * @returns Division object representing @f$ \frac{c}{g} @f$
   */
  template <class Number, class RHSDerived,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<RHSDerived, RealFunction<Number>>>>
  auto
  operator/(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(RealFunction(lhs), rhs);
  }
}
#endif
