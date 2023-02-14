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

#include <mfem.hpp>

#include "Rodin/Geometry/Simplex.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "RangeShape.h"
#include "Exceptions.h"

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
    using Parent = FunctionBase<ScalarFunctionBase<Derived>>;
    public:
      constexpr
      ScalarFunctionBase() = default;

      constexpr
      ScalarFunctionBase(const ScalarFunctionBase& other)
        : Parent(other)
      {}

      constexpr
      ScalarFunctionBase(ScalarFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~ScalarFunctionBase() = default;

      inline
      constexpr
      Scalar getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      constexpr
      RangeShape getRangeShape() const
      {
        return {1, 1};
      }

      virtual ScalarFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup ScalarFunctionSpecializations
   */
  template <class Derived>
  class ScalarFunction<FunctionBase<Derived>>
    : public ScalarFunctionBase<ScalarFunction<FunctionBase<Derived>>>
  {
    using Parent = ScalarFunctionBase<ScalarFunction<FunctionBase<Derived>>>;
    public:
      constexpr
      ScalarFunction(const FunctionBase<Derived>& nested)
        : m_nested(nested)
      {}

      constexpr
      ScalarFunction(const ScalarFunction& other)
        : Parent(other),
          m_nested(other.m_nested)
      {}

      constexpr
      ScalarFunction(ScalarFunction&& other)
        : Parent(std::move(other)),
          m_nested(std::move(other.m_nested))
      {}

      inline
      constexpr
      Scalar getValue(const Geometry::Point& v) const
      {
        return static_cast<Scalar>(m_nested->getValue(v));
      }

      ScalarFunction& traceOf(Geometry::Attribute attrs) override
      {
        ScalarFunctionBase<ScalarFunction<FunctionBase<Derived>>>::traceOf(attrs);
        m_nested->traceOf(attrs);
        return *this;
      }

      ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      FunctionBase<Derived> m_nested;
  };
  template <class Derived>
  ScalarFunction(const FunctionBase<Derived>&) -> ScalarFunction<FunctionBase<Derived>>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a ScalarFunction of arithmetic type `Number`.
   *
   * @tparam Number Arithmetic type
   * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
   */
  template <class Number>
  class ScalarFunction<Number, std::enable_if_t<std::is_arithmetic_v<Number>>>
    : public ScalarFunctionBase<ScalarFunction<Number, std::enable_if_t<std::is_arithmetic_v<Number>>>>
  {
    using Parent =
      ScalarFunctionBase<ScalarFunction<Number, std::enable_if_t<std::is_arithmetic_v<Number>>>>;
    public:
      /**
       * @brief Constructs a ScalarFunction from an arithmetic value.
       * @param[in] x Arithmetic value
       */
      constexpr
      ScalarFunction(Number x)
        : m_x(x)
      {}

      constexpr
      ScalarFunction(const ScalarFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      constexpr
      ScalarFunction(ScalarFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      inline
      constexpr
      Scalar getValue(const Geometry::Point&) const
      {
        return static_cast<Scalar>(m_x);
      }

      ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

      Internal::MFEMFunction build(const Geometry::MeshBase&) const override
      {
        return Internal::MFEMFunction(new mfem::ConstantCoefficient(m_x));
      }

    private:
      const Number m_x;
  };
  template <class Number>
  ScalarFunction(Number)
    -> ScalarFunction<Number, std::enable_if_t<std::is_arithmetic_v<Number>>>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a scalar function given by an arbitrary scalar function.
   */
  template <>
  class ScalarFunction<std::function<Scalar(const Geometry::Point&)>>
    : public ScalarFunctionBase<ScalarFunction<std::function<double(const Geometry::Point&)>>>
  {
    public:
      template <class T>
      constexpr
      ScalarFunction(T&& f)
        : ScalarFunction(std::function<Scalar(const Geometry::Point&)>(std::forward<T>(f)))
      {}

      /**
       * @brief Constructs a ScalarFunction from an std::function.
       */
      ScalarFunction(std::function<Scalar(const Geometry::Point&)> f)
        : m_f(f)
      {}

      ScalarFunction(const ScalarFunction& other)
        : ScalarFunctionBase(other),
          m_f(other.m_f)
      {}

      ScalarFunction(ScalarFunction&& other)
        : ScalarFunctionBase(std::move(other)),
          m_f(std::move(other.m_f))
      {}

      inline
      constexpr
      Scalar getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      const std::function<Scalar(const Geometry::Point&)> m_f;
  };

  ScalarFunction(std::function<Scalar(const Geometry::Point&)>)
    -> ScalarFunction<std::function<Scalar(const Geometry::Point&)>>;

  template <class T>
  ScalarFunction(T)
    -> ScalarFunction<
      std::enable_if_t<std::is_invocable_r_v<Scalar, T, const Geometry::Point&>,
      std::function<Scalar(const Geometry::Point&)>>>;
}

#endif
