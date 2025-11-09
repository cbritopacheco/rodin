/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RealFunction.h
 * @brief Real-valued scalar functions for variational formulations.
 *
 * This file defines RealFunctionBase and RealFunction, representing functions
 * mapping points to real numbers: @f$ f: \Omega \to \mathbb{R} @f$.
 */
#ifndef RODIN_VARIATIONAL_REALFUNCTION_H
#define RODIN_VARIATIONAL_REALFUNCTION_H

#include <memory>
#include <type_traits>

#include "ForwardDecls.h"

#include "ScalarFunction.h"

namespace Rodin::FormLanguage
{
  template <class Derived>
  struct Traits<Variational::RealFunctionBase<Derived>>
  {
    using ScalarType = Real;
    using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup RealFunctionSpecializations RealFunction Template Specializations
   * @brief Template specializations of the RealFunction class.
   * @see RealFunction
   */

  /**
   * @brief Base class for real-valued scalar functions.
   *
   * RealFunctionBase extends ScalarFunctionBase with Real as the scalar type,
   * representing functions @f$ f: \Omega \to \mathbb{R} @f$. This is the most
   * common function type in finite element analysis, used for scalar fields like
   * temperature, pressure, or concentration.
   *
   * @tparam Derived The derived class following CRTP pattern
   *
   * ## Usage Examples
   * ```
   * // Constant real function
   * RealFunction<Real> f(3.14);
   *
   * // Lambda-based function
   * RealFunction<std::function<Real(const Geometry::Point&)>> g(
   *   [](const Geometry::Point& p) { return p.x() * p.y(); }
   * );
   * ```
   *
   * @see ScalarFunctionBase, RealFunction
   */
  template <class Derived>
  class RealFunctionBase : public ScalarFunctionBase<Real, RealFunctionBase<Derived>>
  {
    public:
      /// @brief Type of scalar values (Real)
      using ScalarType = Real;

      /// @brief Parent class type
      using Parent = ScalarFunctionBase<ScalarType, RealFunctionBase<Derived>>;

      /// @brief Import traceOf methods from parent
      using Parent::traceOf;

      /// @brief Import operator() from parent
      using Parent::operator();

      /// @brief Default constructor
      RealFunctionBase() = default;

      /// @brief Copy constructor
      /// @param[in] other Function to copy from
      RealFunctionBase(const RealFunctionBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor
      /// @param[in] other Function to move from
      RealFunctionBase(RealFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /// @brief Virtual destructor
      virtual ~RealFunctionBase() = default;

      /**
       * @brief Evaluates the function at a point.
       *
       * CRTP method delegating to derived class implementation.
       *
       * @param[in] p Point at which to evaluate
       * @returns Real value at the point
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
      virtual RealFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup RealFunctionSpecializations
   */
  template <class NestedDerived>
  class RealFunction<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<FunctionBase<NestedDerived>>
  {
    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<FunctionBase<NestedDerived>>;

      RealFunction(const RealFunctionBase<NestedDerived>& nested)
        : m_nested(nested.copy())
      {}

      RealFunction(const RealFunction& other)
        : Parent(other),
          m_nested(other.m_nested->copy())
      {}

      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_nested(std::move(other.m_nested))
      {}

      constexpr
      decltype(auto) getValue(const Geometry::Point& v) const
      {
        return m_nested->getValue(v);
      }

      RealFunction* copy() const noexcept override
      {
        return new RealFunction(*this);
      }

    private:
      std::unique_ptr<FunctionBase<NestedDerived>> m_nested;
  };

  /**
   * @brief CTAD for RealFunction.
   */
  template <class Derived>
  RealFunction(const RealFunctionBase<Derived>&) -> RealFunction<FunctionBase<Derived>>;

  /**
   * @ingroup RealFunctionSpecializations
   * @brief Represents a constant scalar function with type Real.
   */
  template <>
  class RealFunction<Real> final
    : public RealFunctionBase<RealFunction<Real>>
  {
    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<RealFunction<Real>>;

      /**
       * @brief Constructs a RealFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      RealFunction(const Real& x)
        : m_x(x)
      {}

      RealFunction(const RealFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_x(std::move(other.m_x))
      {}

      constexpr
      Real getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      template <class ... Args>
      RealFunction& traceOf(Args&&... args) noexcept
      {
        return *this;
      }

      RealFunction* copy() const noexcept override
      {
        return new RealFunction(*this);
      }

    private:
      const Real m_x;
  };

  /**
   * @brief CTAD for RealFunction.
   */
  RealFunction(Real) -> RealFunction<Real>;

  template <>
  class RealFunction<Integer> final
    : public RealFunctionBase<RealFunction<Integer>>
  {
    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<RealFunction<Integer>>;

      /**
       * @brief Constructs a RealFunction from an integer value.
       * @param[in] x Constant integer value
       */
      RealFunction(Integer x)
        : m_x(x)
      {}

      RealFunction(const RealFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_x(std::move(other.m_x))
      {}

      constexpr
      Real getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      template <class ... Args>
      RealFunction& traceOf(Args&&... args) noexcept
      {
        return *this;
      }

      RealFunction* copy() const noexcept override
      {
        return new RealFunction(*this);
      }

    private:
      const Integer m_x;
  };

  RealFunction(Integer) -> RealFunction<Integer>;

  /**
   * @ingroup RealFunctionSpecializations
   * @brief Represents a scalar function given by an arbitrary scalar function.
   */
  template <class F>
  class RealFunction<F> final : public RealFunctionBase<RealFunction<F>>
  {
    static_assert(std::is_invocable_r_v<Real, F, const Geometry::Point&>);

    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<RealFunction<F>>;

      RealFunction(F f)
        : m_f(f)
      {}

      RealFunction(const RealFunction& other)
        : Parent(other),
          m_f(other.m_f)
      {}

      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f))
      {}

      constexpr
      Real getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      template <class ... Args>
      RealFunction& traceOf(Args&&... args) noexcept
      {
        return *this;
      }

      RealFunction* copy() const noexcept override
      {
        return new RealFunction(*this);
      }

    private:
      const F m_f;
  };

  /**
   * @brief CTAD for RealFunction.
   */
  template <class F, typename = std::enable_if_t<std::is_invocable_v<F, const Geometry::Point&>>>
  RealFunction(F) -> RealFunction<F>;
}

#endif
