/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_AND_H
#define RODIN_VARIATIONAL_AND_H

/**
 * @file
 * @brief Logical AND operator for boolean functions.
 */

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ANDSpecializations AND Template Specializations
   * @brief Template specializations of the AND class.
   * @see AND
   */

  /**
   * @ingroup ANDSpecializations
   * @brief Logical AND operator between two boolean functions.
   *
   * Returns a boolean function that evaluates to the logical conjunction:
   * @f[
   *    \text{AND}(f, g)(x) = f(x) \land g(x)
   * @f]
   *
   * Common applications:
   * - Region intersection: points that satisfy multiple conditions
   * - Composite boundary conditions
   * - Multi-criteria selection
   *
   * @tparam LHSDerived Left-hand side boolean function type
   * @tparam RHSDerived Right-hand side boolean function type
   *
   * @see operator&&
   */
  template <class LHSDerived, class RHSDerived>
  class AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = BooleanFunctionBase<LHSDerived>;

      using RHSType = BooleanFunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<AND<LHSType, RHSType>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs the logical AND of two boolean functions.
       * @param lhs Left-hand side boolean function
       * @param rhs Right-hand side boolean function
       */
      AND(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy
       */
      AND(const AND& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move
       */
      AND(AND&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Gets the left-hand side function.
       * @returns Reference to left-hand side function
       */
      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      /**
       * @brief Gets the right-hand side function.
       * @returns Reference to right-hand side function
       */
      const auto& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      /**
       * @brief Evaluates the logical AND at a point.
       * @param p Point at which to evaluate
       * @returns Logical conjunction of operands at p
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) && getRHS().getValue(p);
      }

      /**
       * @brief Copies the AND object.
       * @returns Pointer to copied object
       */
      AND* copy() const noexcept final override
      {
        return new AND(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief Deduction guide for AND.
   */
  template <class LHSDerived, class RHSDerived>
  AND(const BooleanFunctionBase<LHSDerived>&, const BooleanFunctionBase<RHSDerived>&)
    -> AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>;

  /**
   * @brief Logical AND operator for boolean functions.
   * @param lhs Left-hand side boolean function
   * @param rhs Right-hand side boolean function
   * @returns AND object representing @f$ lhs \land rhs @f$
   */
  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator&&(const BooleanFunctionBase<LHSDerived>& lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return AND(lhs, rhs);
  }

  /**
   * @brief Logical AND with boolean constant on left.
   * @param lhs Boolean constant
   * @param rhs Boolean function
   * @returns AND object
   */
  template <class RHSDerived>
  constexpr
  auto
  operator&&(Boolean lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return AND(BooleanFunction(lhs), rhs);
  }

  /**
   * @brief Logical AND with boolean constant on right.
   * @param lhs Boolean function
   * @param rhs Boolean constant
   * @returns AND object
   */
  template <class LHSDerived>
  constexpr
  auto
  operator&&(const BooleanFunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return AND(lhs, BooleanFunction(rhs));
  }
}

#endif


