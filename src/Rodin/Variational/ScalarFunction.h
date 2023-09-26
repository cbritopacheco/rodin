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

#include "Function.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ScalarFunctionSpecializations ScalarFunction Template Specializations
   * @brief Template specializations of the ScalarFunction class.
   * @see ScalarFunction
   */

  template <class Derived>
  class ScalarFunctionBase : public FunctionBase<ScalarFunctionBase<Derived>>
  {
    public:
      using Parent = FunctionBase<ScalarFunctionBase<Derived>>;

      ScalarFunctionBase() = default;

      ScalarFunctionBase(const ScalarFunctionBase& other)
        : Parent(other)
      {}

      ScalarFunctionBase(ScalarFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~ScalarFunctionBase() = default;

      inline
      constexpr
      ScalarFunctionBase& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      ScalarFunctionBase& traceOf(const std::set<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

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

      inline
      constexpr
      void getValue(Math::Vector&, const Geometry::Point&) const = delete;

      inline
      constexpr
      void getValue(Math::Matrix&, const Geometry::Point&) const = delete;

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { 1, 1 };
      }

      virtual ScalarFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };

  /**
   * @ingroup ScalarFunctionSpecializations
   */
  template <class NestedDerived>
  class ScalarFunction<FunctionBase<NestedDerived>> final
    : public ScalarFunctionBase<ScalarFunction<NestedDerived>>
  {
    public:
      using Parent = ScalarFunctionBase<ScalarFunction<NestedDerived>>;
      using NestedRangeType = typename FormLanguage::Traits<FunctionBase<NestedDerived>>::RangeType;

      static_assert(std::is_same_v<NestedRangeType, Scalar>);

      ScalarFunction(const FunctionBase<NestedDerived>& nested)
        : m_nested(nested.copy())
      {}

      ScalarFunction(const ScalarFunction& other)
        : Parent(other),
          m_nested(other.m_nested->copy())
      {}

      ScalarFunction(ScalarFunction&& other)
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
      ScalarFunction& traceOf(Geometry::Attribute attrs)
      {
        m_nested->traceOf(attrs);
        return *this;
      }

      inline ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      std::unique_ptr<FunctionBase<NestedDerived>> m_nested;
  };

  template <class Derived>
  ScalarFunction(const FunctionBase<Derived>&) -> ScalarFunction<FunctionBase<Derived>>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a constant scalar function with type Scalar.
   */
  template <>
  class ScalarFunction<Scalar> final
    : public ScalarFunctionBase<ScalarFunction<Scalar>>
  {
    public:
      using Parent = ScalarFunctionBase<ScalarFunction<Scalar>>;

      /**
       * @brief Constructs a ScalarFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      ScalarFunction(Scalar x)
        : m_x(x)
      {}

      ScalarFunction(const ScalarFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      ScalarFunction(ScalarFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      inline
      constexpr
      ScalarFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Scalar getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      Scalar getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      inline ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      const Scalar m_x;
  };

  ScalarFunction(Scalar) -> ScalarFunction<Scalar>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a constant scalar function with type Integer.
   */
  template <>
  class ScalarFunction<Integer> final
    : public ScalarFunctionBase<ScalarFunction<Integer>>
  {
    public:
      using Parent = ScalarFunctionBase<ScalarFunction<Integer>>;

      /**
       * @brief Constructs a ScalarFunction from an integer value.
       * @param[in] x Constant integer value
       */
      ScalarFunction(Integer x)
        : m_x(x)
      {}

      ScalarFunction(const ScalarFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      ScalarFunction(ScalarFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      inline
      constexpr
      ScalarFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Scalar getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      Scalar getValue(const Geometry::Point&) const
      {
        return static_cast<Scalar>(m_x);
      }

      inline ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      const Integer m_x;
  };

  ScalarFunction(Integer) -> ScalarFunction<Integer>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a scalar function given by an arbitrary scalar function.
   */
  template <class F>
  class ScalarFunction<F> final : public ScalarFunctionBase<ScalarFunction<F>>
  {
    static_assert(std::is_invocable_r_v<Scalar, F, const Geometry::Point&>);

    public:
      using Parent = ScalarFunctionBase<ScalarFunction<F>>;

      ScalarFunction(F f)
        : m_f(f)
      {}

      ScalarFunction(const ScalarFunction& other)
        : Parent(other),
          m_f(other.m_f)
      {}

      ScalarFunction(ScalarFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f))
      {}

      inline
      constexpr
      ScalarFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Scalar getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      inline ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      const F m_f;
  };

  template <class F, typename = std::enable_if_t<std::is_invocable_r_v<Scalar, F, const Geometry::Point&>>>
  ScalarFunction(F) -> ScalarFunction<F>;

  namespace P
  {
    static ScalarFunction x([](const Geometry::Point& p) { return p.x(); });
    static ScalarFunction y([](const Geometry::Point& p) { return p.y(); });
    static ScalarFunction z([](const Geometry::Point& p) { return p.z(); });
  }
}

#endif
