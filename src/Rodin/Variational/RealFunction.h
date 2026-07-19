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
  /// @brief Traits for the CRTP base of real-valued functions.
  template <class Derived>
  struct Traits<Variational::RealFunctionBase<Derived>>
  {
    /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Derived CRTP function type.
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
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Evaluates the function at an integration point.
       * @param[in] ip Integration point at which to evaluate
       * @returns Real value at the integration point
       */
      constexpr
      auto getValue(const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, const IntegrationPoint& q) { f.getValue(q); })
          return static_cast<const Derived&>(*this).getValue(ip);
        else
          return static_cast<const Derived&>(*this).getValue(ip.getPoint());
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
       * @brief Returns the polynomial order on a mesh entity.
       * @param[in] geom Entity whose local polynomial order is requested
       * @returns Polynomial order when known
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(geom);
      }

      /**
       * @brief Creates a polymorphic copy of the function.
       * @returns Pointer to newly allocated copy
       */
      virtual RealFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup RealFunctionSpecializations
   * @brief Real-valued function that wraps another real function expression.
   */
  template <class NestedDerived>
  class RealFunction<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<FunctionBase<NestedDerived>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Parent class type.
      using Parent = RealFunctionBase<FunctionBase<NestedDerived>>;

      /**
       * @brief Constructs a real function wrapper around an existing expression.
       * @param[in] nested Real-valued expression to clone
       */
      RealFunction(const RealFunctionBase<NestedDerived>& nested)
        : m_nested(nested.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Function to copy
       */
      RealFunction(const RealFunction& other)
        : Parent(other),
          m_nested(other.m_nested->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Function to move from
       */
      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_nested(std::move(other.m_nested))
      {}

      /**
       * @brief Evaluates the wrapped function at a point.
       * @param[in] v Point at which to evaluate
       * @returns Wrapped function value
       */
      constexpr
      auto getValue(const Geometry::Point& v) const
      {
        return m_nested->getValue(v);
      }

      /**
       * @brief Evaluates the wrapped function at an integration point.
       * @param[in] ip Integration point at which to evaluate
       * @returns Wrapped function value
       */
      constexpr
      auto getValue(const IntegrationPoint& ip) const
      {
        return m_nested->getValue(ip);
      }

      /**
       * @brief Returns the wrapped function order on a mesh entity.
       * @param[in] geom Entity whose local polynomial order is requested
       * @returns Polynomial order when known
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return m_nested->getOrder(geom);
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
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Parent class type.
      using Parent = RealFunctionBase<RealFunction<Real>>;

      /**
       * @brief Constructs a RealFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      RealFunction(const Real& x)
        : m_x(x)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Function to copy
       */
      RealFunction(const RealFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Function to move from
       */
      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_x(std::move(other.m_x))
      {}

      /**
       * @brief Evaluates the constant function at a point.
       * @returns Constant real value
       */
      constexpr
      Real getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      /**
       * @brief Leaves constant functions unchanged on trace domains.
       * @returns Reference to this function
       */
      template <class ... Args>
      RealFunction& traceOf(Args&&... args) noexcept
      {
        return *this;
      }

      /**
       * @brief Returns the polynomial order of a constant function.
       * @returns Zero polynomial order
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return size_t(0);
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

  /**
   * @ingroup RealFunctionSpecializations
   * @brief Represents a constant scalar function initialized from an Integer.
   */
  template <>
  class RealFunction<Integer> final
    : public RealFunctionBase<RealFunction<Integer>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Parent class type.
      using Parent = RealFunctionBase<RealFunction<Integer>>;

      /**
       * @brief Constructs a RealFunction from an integer value.
       * @param[in] x Constant integer value
       */
      RealFunction(Integer x)
        : m_x(x)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Function to copy
       */
      RealFunction(const RealFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Function to move from
       */
      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_x(std::move(other.m_x))
      {}

      /**
       * @brief Evaluates the constant function at a point.
       * @returns Constant value converted to Real
       */
      constexpr
      Real getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      /**
       * @brief Leaves constant functions unchanged on trace domains.
       * @returns Reference to this function
       */
      template <class ... Args>
      RealFunction& traceOf(Args&&... args) noexcept
      {
        return *this;
      }

      /**
       * @brief Returns the polynomial order of a constant function.
       * @returns Zero polynomial order
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return size_t(0);
      }

      RealFunction* copy() const noexcept override
      {
        return new RealFunction(*this);
      }

    private:
      const Integer m_x;
  };

  /// @brief CTAD for integer constants.
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
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Parent class type.
      using Parent = RealFunctionBase<RealFunction<F>>;

      /**
       * @brief Constructs a real function from a callable.
       * @param[in] f Callable returning a Real from a Geometry::Point
       */
      RealFunction(F f)
        : m_f(f)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Function to copy
       */
      RealFunction(const RealFunction& other)
        : Parent(other),
          m_f(other.m_f)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Function to move from
       */
      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f))
      {}

      /**
       * @brief Evaluates the callable at a point.
       * @param[in] v Point at which to evaluate
       * @returns Callable value
       */
      constexpr
      Real getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      /**
       * @brief Evaluates the callable at an integration point.
       * @param[in] ip Integration point at which to evaluate
       * @returns Callable value
       */
      constexpr
      Real getValue(const IntegrationPoint& ip) const
      {
        if constexpr (std::is_invocable_r_v<Real, F, const IntegrationPoint&>)
          return m_f(ip);
        else
          return m_f(ip.getPoint());
      }

      /**
       * @brief Leaves callable functions unchanged on trace domains.
       * @returns Reference to this function
       */
      template <class ... Args>
      RealFunction& traceOf(Args&&... args) noexcept
      {
        return *this;
      }

      /**
       * @brief Returns the polynomial order of the callable.
       * @returns std::nullopt because arbitrary callables have unknown order
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
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
