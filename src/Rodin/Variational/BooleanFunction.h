/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BooleanFunction.h
 * @brief Boolean-valued functions for logical operations in variational formulations.
 *
 * This file defines BooleanFunctionBase and BooleanFunction for representing
 * functions mapping points to boolean values: @f$ f: \Omega \to \{0, 1\} @f$.
 * These are commonly used for characteristic functions, indicator functions,
 * and logical conditions in variational problems.
 */
#ifndef RODIN_VARIATIONAL_BOOLEANFUNCTION_H
#define RODIN_VARIATIONAL_BOOLEANFUNCTION_H

#include "ForwardDecls.h"
#include "Function.h"

namespace Rodin::Variational
{
  /**
   * @defgroup BooleanFunctionSpecializations BooleanFunction Template Specializations
   * @brief Template specializations of the BooleanFunction class.
   * @see BooleanFunction
   */

  /**
   * @brief Base class for boolean-valued functions.
   *
   * BooleanFunctionBase represents functions mapping points to boolean values:
   * @f$ f: \Omega \to \{0, 1\} @f$ or equivalently @f$ f: \Omega \to \{\text{true}, \text{false}\} @f$.
   *
   * These functions are commonly used for:
   * - **Indicator functions**: @f$ \chi_A(x) = \begin{cases} 1 & x \in A \\ 0 & x \notin A \end{cases} @f$
   * - **Logical conditions**: Specifying regions for boundary conditions or material properties
   * - **Comparison operators**: Results of @f$ f(x) > g(x) @f$, @f$ f(x) = g(x) @f$, etc.
   *
   * @tparam Derived The derived class following CRTP pattern
   *
   * ## Usage Examples
   * ```
   * // Constant boolean function
   * BooleanFunction<Boolean> isActive(true);
   *
   * // Indicator function for a region
   * auto indicator = (x > 0) && (y < 1);  // Returns BooleanFunction
   * ```
   *
   * @see FunctionBase, AND, OR, EQ, GT, LT
   */
  template <class Derived>
  class BooleanFunctionBase
    : public FunctionBase<BooleanFunctionBase<Derived>>
  {
    public:
      /// @brief Parent class type
      using Parent = FunctionBase<BooleanFunctionBase<Derived>>;
      
      /// @brief Import operator() from parent
      using Parent::operator();

      /// @brief Default constructor
      BooleanFunctionBase() = default;

      /// @brief Copy constructor
      /// @param[in] other Function to copy from
      BooleanFunctionBase(const BooleanFunctionBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor
      /// @param[in] other Function to move from
      BooleanFunctionBase(BooleanFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /// @brief Virtual destructor
      virtual ~BooleanFunctionBase() = default;

      /**
       * @brief Evaluates the boolean function at a point.
       *
       * CRTP method that delegates to the derived class's implementation.
       *
       * @param[in] p Point at which to evaluate
       * @returns Boolean value at the given point
       * @note CRTP function to be overridden in Derived class.
       */
      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Sets the trace domain for the function.
       *
       * @tparam Args Variadic template for trace domain specification
       * @param[in] args Arguments specifying the trace domain
       * @returns Reference to derived object (for method chaining)
       */
      template <class ... Args>
      constexpr
      Derived& traceOf(const Args& ... args)
      {
        return static_cast<Derived&>(*this).traceOf(args...);
      }

      /**
       * @brief Creates a polymorphic copy of the function.
       * @returns Pointer to newly allocated copy
       */
      virtual BooleanFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup BooleanFunctionSpecializations
   */
  template <>
  class BooleanFunction<Boolean> final
    : public BooleanFunctionBase<BooleanFunction<Boolean>>
  {
    public:
      using Parent = BooleanFunctionBase<BooleanFunction<Boolean>>;

      BooleanFunction(Boolean v)
        : m_v(v)
      {}

      BooleanFunction(const BooleanFunction& other)
        : Parent(other),
          m_v(other.m_v)
      {}

      BooleanFunction(BooleanFunction&& other)
        : Parent(std::move(other)),
          m_v(other.m_v)
      {}

      constexpr
      Boolean getValue(const Geometry::Point&) const
      {
        return m_v;
      }

      template <class ... Args>
      constexpr
      BooleanFunction& traceOf(const Args& ... args)
      {
        return *this;
      }

      BooleanFunction* copy() const noexcept final override
      {
        return new BooleanFunction(*this);
      }

    private:
      const Boolean m_v;
  };

  BooleanFunction(Boolean) -> BooleanFunction<Boolean>;
}

#endif

