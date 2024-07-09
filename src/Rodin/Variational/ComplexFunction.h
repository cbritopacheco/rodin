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
   * @defgroup ComplexFunctionSpecializations ComplexFunction Template Specializations
   * @brief Template specializations of the ComplexFunction class.
   * @see ComplexFunction
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
      ComplexFunction(const Integer& x)
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
      const Integer& getValue() const
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

  template <>
  class ComplexFunction<Real> final
    : public ComplexFunctionBase<ComplexFunction<Real>>
  {
    public:
      using Parent = ComplexFunctionBase<ComplexFunction<Real>>;

      /**
       * @brief Constructs a ComplexFunction from an Real value.
       * @param[in] x Constant Real value
       */
      ComplexFunction(const Real& x)
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
      const Real& getValue() const
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
      const Real m_x;
  };

  ComplexFunction(Real) -> ComplexFunction<Real>;

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
      ComplexFunction(const Complex& x)
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
      const Complex& getValue() const
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
  ComplexFunction(const Complex&) -> ComplexFunction<Complex>;

  /**
   * @ingroup ComplexFunctionSpecializations
   */
  template <class NestedDerived>
  class ComplexFunction<FunctionBase<NestedDerived>> final
    : public ComplexFunctionBase<ComplexFunction<NestedDerived>>
  {
    public:
      using Parent = ComplexFunctionBase<ComplexFunction<NestedDerived>>;

      using ScalarType = Complex;

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
      ScalarType getValue(const Geometry::Point& v) const
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
   */
  template <class RealNestedDerived, class ImagNestedDerived>
  class ComplexFunction<FunctionBase<RealNestedDerived>, FunctionBase<ImagNestedDerived>> final
    : public ComplexFunctionBase<ComplexFunction<FunctionBase<RealNestedDerived>, FunctionBase<ImagNestedDerived>>>
  {
    public:
      using RealFunctionType = FunctionBase<RealNestedDerived>;

      using ImagFunctionType = FunctionBase<ImagNestedDerived>;

      using Parent =
        ComplexFunctionBase<ComplexFunction<RealFunctionType, ImagFunctionType>>;

      using RealFunctionRangeType = typename FormLanguage::Traits<RealFunctionType>::RangeType;

      using ImagFunctionRangeType = typename FormLanguage::Traits<ImagFunctionType>::RangeType;

      static_assert(std::is_same_v<RealFunctionType, Real>);

      static_assert(std::is_same_v<ImagFunctionType, Real>);

      ComplexFunction(const RealFunctionType& re, const ImagFunctionType& imag)
        : m_re(re.copy()), m_imag(imag.copy())
      {}

      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_re(m_re->copy()), m_imag(m_imag->copy())
      {}

      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_re(std::move(m_re)), m_imag(std::move(m_imag))
      {}

      inline
      constexpr
      Complex getValue(const Geometry::Point& p) const
      {
        return { m_re->getValue(p), m_imag->getValue(p) };
      }

      inline
      constexpr
      ComplexFunction& traceOf(Geometry::Attribute attrs)
      {
        m_re->traceOf(attrs);
        m_imag->traceOf(attrs);
        return *this;
      }

      inline ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      std::unique_ptr<RealFunctionType> m_re;
      std::unique_ptr<RealFunctionType> m_imag;
  };

  /**
   * @brief CTAD for ComplexFunction.
   */
  template <class RealNestedDerived, class ImagNestedDerived>
  ComplexFunction(const FunctionBase<RealNestedDerived>&, const FunctionBase<ImagNestedDerived>&)
    -> ComplexFunction<FunctionBase<RealNestedDerived>, FunctionBase<ImagNestedDerived>>;

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

      ComplexFunction(const F& f)
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
  ComplexFunction(const F&) -> ComplexFunction<F>;

  template <class FReal, class FImag>
  class ComplexFunction<FReal, FImag> final
    : public ComplexFunctionBase<ComplexFunction<FReal, FImag>>
  {
    static_assert(std::is_invocable_v<FReal, const Geometry::Point&>);
    static_assert(std::is_invocable_v<FImag, const Geometry::Point&>);

    public:
      using Parent = ComplexFunctionBase<ComplexFunction<FReal, FImag>>;

      ComplexFunction(const FReal& re, const FImag& imag)
        : m_re(re), m_imag(imag)
      {}

      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_re(other.m_re), m_imag(other.m_imag)
      {}

      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_re(std::move(other.m_re)), m_imag(std::move(other.m_imag))
      {}

      inline
      constexpr
      ComplexFunction& traceOf(Geometry::Attribute attr)
      {
        m_re->traceOf(attr);
        m_imag->traceOf(attr);
        return *this;
      }

      inline
      constexpr
      ComplexFunction& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_re->traceOf(attrs);
        m_imag->traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      Complex getValue(const Geometry::Point& p) const
      {
        return { m_re->getValue(p), m_imag->getValue(p) };
      }

      inline ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      const FReal m_re;
      const FImag m_imag;
  };

  /**
   * @brief CTAD for ComplexFunction.
   */
  template <class FReal, class FImag, typename =
    std::enable_if_t<
      std::is_invocable_v<FReal, const Geometry::Point&> && std::is_invocable_v<FImag, const Geometry::Point&>>>
  ComplexFunction(const FReal&, const FImag&) -> ComplexFunction<FReal, FImag>;
}

#endif

