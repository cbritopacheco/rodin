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

#include "Function.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Derived>
  struct Traits<Variational::ScalarFunctionBase<Scalar, Derived>>
  {
    using ScalarType = Scalar;
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

  template <class Scalar, class Derived>
  class ScalarFunctionBase
    : public FunctionBase<ScalarFunctionBase<Scalar, Derived>>
  {
    public:
      using ScalarType = Scalar;

      using Parent = FunctionBase<ScalarFunctionBase<ScalarType, Derived>>;

      using Parent::traceOf;

      ScalarFunctionBase() = default;

      ScalarFunctionBase(const ScalarFunctionBase& other)
        : Parent(other)
      {}

      ScalarFunctionBase(ScalarFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~ScalarFunctionBase() = default;

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
      void getValue(ScalarType& res, const Geometry::Point& p) const
      {
        res = this->getValue(p);
      }

      virtual ScalarFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup ScalarFunctionSpecializations
   */
  template <class Scalar, class NestedDerived>
  class ScalarFunction<ScalarFunctionBase<Scalar, NestedDerived>> final
    : public ScalarFunctionBase<Scalar, ScalarFunction<ScalarFunctionBase<Scalar, NestedDerived>>>
  {
    public:
      using ScalarType = Scalar;

      using Parent =
        ScalarFunctionBase<ScalarType, ScalarFunction<ScalarFunctionBase<ScalarType, NestedDerived>>>;

      ScalarFunction(const ScalarFunctionBase<ScalarType, NestedDerived>& nested)
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
      std::unique_ptr<ScalarFunctionBase<ScalarType, NestedDerived>> m_nested;
  };

  /**
   * @brief CTAD for ScalarFunction.
   */
  template <class Scalar, class NestedDerived>
  ScalarFunction(const ScalarFunctionBase<Scalar, NestedDerived>&)
    -> ScalarFunction<ScalarFunctionBase<Scalar, NestedDerived>>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a constant scalar function with type Real.
   */
  template <>
  class ScalarFunction<Real> final
    : public ScalarFunctionBase<Real, ScalarFunction<Real>>
  {
    public:
      using ScalarType = Real;

      using Parent = ScalarFunctionBase<Real, ScalarFunction<Real>>;

      /**
       * @brief Constructs a ScalarFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      ScalarFunction(const ScalarType& x)
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
      ScalarType getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      ScalarType getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      inline ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      const ScalarType m_x;
  };

  /**
   * @brief CTAD for ScalarFunction.
   */
  ScalarFunction(const Real&) -> ScalarFunction<Real>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a constant scalar function with type Complex.
   */
  template <>
  class ScalarFunction<Complex> final
    : public ScalarFunctionBase<Complex, ScalarFunction<Complex>>
  {
    public:
      using ScalarType = Complex;

      using Parent = ScalarFunctionBase<Complex, ScalarFunction<Complex>>;

      /**
       * @brief Constructs a ScalarFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      ScalarFunction(const ScalarType& x)
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
      ScalarType getValue() const
      {
        return m_x;
      }

      inline
      constexpr
      ScalarType getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      inline ScalarFunction* copy() const noexcept override
      {
        return new ScalarFunction(*this);
      }

    private:
      const ScalarType m_x;
  };

  /**
   * @brief CTAD for ScalarFunction.
   */
  ScalarFunction(const Complex&) -> ScalarFunction<Complex>;

  /**
   * @ingroup ScalarFunctionSpecializations
   * @brief Represents a scalar function given by an arbitrary scalar function.
   */
  template <class F>
  class ScalarFunction<F> final
    : public ScalarFunctionBase<std::invoke_result_t<F, const Geometry::Point&>, ScalarFunction<F>>
  {

    public:
      using ScalarType = std::invoke_result_t<F, const Geometry::Point&>;

      using Parent = ScalarFunctionBase<ScalarType, ScalarFunction<F>>;

      static_assert(std::is_invocable_v<F, const Geometry::Point&>);

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
      ScalarType getValue(const Geometry::Point& v) const
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

  /**
   * @brief CTAD for ScalarFunction.
   */
  template <class F, typename = std::enable_if_t<std::is_invocable_v<F, const Geometry::Point&>>>
  ScalarFunction(F) -> ScalarFunction<F>;
}

#endif

