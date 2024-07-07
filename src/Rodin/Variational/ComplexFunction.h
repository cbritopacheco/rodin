/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COMPLEXFUNCTION_H
#define RODIN_VARIATIONAL_COMPLEXFUNCTION_H

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
  class ComplexFunctionBase : public ScalarFunctionBase<Complex, ComplexFunctionBase<Derived>>
  {
    public:
      using ScalarType = Complex;

      using Parent = ScalarFunctionBase<ScalarType, ComplexFunctionBase<Derived>>;

      using Parent::traceOf;

      ComplexFunctionBase() = default;

      ComplexFunctionBase(const ComplexFunctionBase& other)
        : Parent(other)
      {}

      ComplexFunctionBase(ComplexFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~ComplexFunctionBase() = default;

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

      virtual ComplexFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup ComplexFunctionSpecializations
   */
  template <class NestedDerived>
  class ComplexFunction<FunctionBase<NestedDerived>> final
    : public ComplexFunctionBase<ComplexFunction<NestedDerived>>
  {
    public:
      using Parent = ComplexFunctionBase<ComplexFunction<NestedDerived>>;

      using NestedRangeType = typename FormLanguage::Traits<FunctionBase<NestedDerived>>::RangeType;

      static_assert(std::is_same_v<NestedRangeType, Complex>);

      ComplexFunction(const FunctionBase<NestedDerived>& nested)
        : m_nested(nested.copy())
      {}

      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_nested(other.m_nested->copy())
      {}

      ComplexFunction(ComplexFunction&& other)
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
      ComplexFunction& traceOf(Geometry::Attribute attrs)
      {
        m_nested->traceOf(attrs);
        return *this;
      }

      inline ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      std::unique_ptr<FunctionBase<NestedDerived>> m_nested;
  };

  /**
   * @brief CTAD for ComplexFunction.
   */
  template <class Derived>
  ComplexFunction(const FunctionBase<Derived>&) -> ComplexFunction<FunctionBase<Derived>>;

  /**
   * @ingroup ComplexFunctionSpecializations
   * @brief Represents a constant scalar function with type Complex.
   */
  template <>
  class ComplexFunction<Complex> final
    : public ComplexFunctionBase<ComplexFunction<Complex>>
  {
    public:
      using Parent = ComplexFunctionBase<ComplexFunction<Complex>>;

      /**
       * @brief Constructs a ComplexFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      ComplexFunction(Complex x)
        : m_x(x)
      {}

      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      inline
      constexpr
      ComplexFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Complex getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      Complex getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      inline ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      const Complex m_x;
  };

  /**
   * @brief CTAD for ComplexFunction.
   */
  ComplexFunction(Complex) -> ComplexFunction<Complex>;

  template <>
  class ComplexFunction<Integer> final
    : public ComplexFunctionBase<ComplexFunction<Integer>>
  {
    public:
      using Parent = ComplexFunctionBase<ComplexFunction<Integer>>;

      /**
       * @brief Constructs a ComplexFunction from an integer value.
       * @param[in] x Constant integer value
       */
      ComplexFunction(Integer x)
        : m_x(x)
      {}

      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      inline
      constexpr
      ComplexFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Complex getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      Complex getValue(const Geometry::Point&) const
      {
        return static_cast<Complex>(m_x);
      }

      inline ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      const Integer m_x;
  };

  ComplexFunction(Integer) -> ComplexFunction<Integer>;

  /**
   * @ingroup ComplexFunctionSpecializations
   * @brief Represents a scalar function given by an arbitrary scalar function.
   */
  template <class F>
  class ComplexFunction<F> final : public ComplexFunctionBase<ComplexFunction<F>>
  {
    static_assert(std::is_invocable_r_v<Complex, F, const Geometry::Point&>);

    public:
      using Parent = ComplexFunctionBase<ComplexFunction<F>>;

      ComplexFunction(F f)
        : m_f(f)
      {}

      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_f(other.m_f)
      {}

      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f))
      {}

      inline
      constexpr
      ComplexFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      Complex getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      inline ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      const F m_f;
  };

  /**
   * @brief CTAD for ComplexFunction.
   */
  template <class F, typename =
    std::enable_if_t<std::is_invocable_r_v<Complex, F, const Geometry::Point&>>>
  ComplexFunction(F) -> ComplexFunction<F>;
}

#endif

