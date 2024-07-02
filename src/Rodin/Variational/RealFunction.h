/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARFUNCTION_H
#define RODIN_VARIATIONAL_SCALARFUNCTION_H

#include <map>
#include <set>
#include <memory>
#include <optional>
#include <type_traits>

#include "Rodin/Geometry/Polytope.h"

#include "ForwardDecls.h"

#include "RangeShape.h"

#include "ScalarFunction.h"


namespace Rodin::Variational
{
  /**
   * @defgroup RealFunctionSpecializations RealFunction Template Specializations
   * @brief Template specializations of the RealFunction class.
   * @see RealFunction
   */

  template <class Derived>
  class RealFunctionBase : public ScalarFunctionBase<Real, RealFunctionBase<Derived>>
  {
    public:
      using ScalarType = Real;

      using Parent = ScalarFunctionBase<ScalarType, RealFunctionBase<Derived>>;

      using Parent::traceOf;

      RealFunctionBase() = default;

      RealFunctionBase(const RealFunctionBase& other)
        : Parent(other)
      {}

      RealFunctionBase(RealFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~RealFunctionBase() = default;

      inline
      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      virtual RealFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup RealFunctionSpecializations
   */
  template <class NestedDerived>
  class RealFunction<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<RealFunction<NestedDerived>>
  {
    public:
      using Parent = RealFunctionBase<RealFunction<NestedDerived>>;

      using NestedRangeType = typename FormLanguage::Traits<FunctionBase<NestedDerived>>::RangeType;

      static_assert(std::is_same_v<NestedRangeType, Real>);

      RealFunction(const FunctionBase<NestedDerived>& nested)
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

      inline
      constexpr
      auto getValue(const Geometry::Point& v) const
      {
        return m_nested->getValue(v);
      }

      inline
      constexpr
      RealFunction& traceOf(Geometry::Attribute attrs)
      {
        m_nested->traceOf(attrs);
        return *this;
      }

      inline RealFunction* copy() const noexcept override
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
  RealFunction(const FunctionBase<Derived>&) -> RealFunction<FunctionBase<Derived>>;

  /**
   * @ingroup RealFunctionSpecializations
   * @brief Represents a constant scalar function with type Real.
   */
  template <>
  class RealFunction<Real> final
    : public RealFunctionBase<RealFunction<Real>>
  {
    public:
      using Parent = RealFunctionBase<RealFunction<Real>>;

      /**
       * @brief Constructs a RealFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      RealFunction(Real x)
        : m_x(x)
      {}

      RealFunction(const RealFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      RealFunction(RealFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      inline
      constexpr
      RealFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Real getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      Real getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      inline RealFunction* copy() const noexcept override
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
          m_x(other.m_x)
      {}

      inline
      constexpr
      RealFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Real getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      Real getValue(const Geometry::Point&) const
      {
        return static_cast<Real>(m_x);
      }

      inline RealFunction* copy() const noexcept override
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

      inline
      constexpr
      RealFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Real getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      inline RealFunction* copy() const noexcept override
      {
        return new RealFunction(*this);
      }

    private:
      const F m_f;
  };

  /**
   * @brief CTAD for RealFunction.
   */
  template <class F, typename =
    std::enable_if_t<std::is_invocable_r_v<Real, F, const Geometry::Point&>>>
  RealFunction(F) -> RealFunction<F>;
}

#endif
