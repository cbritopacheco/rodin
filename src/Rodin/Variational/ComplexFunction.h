/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COMPLEXFUNCTION_H
#define RODIN_VARIATIONAL_COMPLEXFUNCTION_H

/**
 * @file ComplexFunction.h
 * @brief Complex-valued scalar functions.
 */

#include <memory>
#include <type_traits>

#include "ForwardDecls.h"

#include "ScalarFunction.h"

/// @cond RODIN_DOXYGEN_INTERNAL
namespace Rodin::Variational
{
  /**
   * @defgroup ComplexFunctionSpecializations ComplexFunction Template Specializations
   * @brief Template specializations of the ComplexFunction class.
   * @see ComplexFunction
   */

  /**
   * @brief Base class for complex-valued scalar functions.
   *
   * Represents functions @f$ f: \Omega \to \mathbb{C} @f$ where @f$ \mathbb{C} @f$ 
   * denotes the field of complex numbers. Complex functions can be:
   * - Constructed from real and imaginary parts
   * - Decomposed using Re() and Im() operations
   * - Conjugated using Conjugate() operation
   *
   * @tparam Derived The derived class type (CRTP pattern)
   *
   * @see ComplexFunction, Re, Im, Conjugate
   */
  template <class Derived>
  class ComplexFunctionBase : public ScalarFunctionBase<Complex, ComplexFunctionBase<Derived>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Complex;

      /// @brief Parent class type.
      using Parent = ScalarFunctionBase<ScalarType, ComplexFunctionBase<Derived>>;

      using Parent::traceOf;

      using Parent::operator();

      /// @brief Default constructor.
      ComplexFunctionBase() = default;

      /// @brief Copy constructor.
      ComplexFunctionBase(const ComplexFunctionBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      ComplexFunctionBase(ComplexFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /// @brief Destructor.
      virtual ~ComplexFunctionBase() = default;

      /// @brief Evaluates the complex function at a geometry point.
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /// @brief Evaluates the complex function at an integration point.
      constexpr
      auto getValue(const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, const IntegrationPoint& q) { f.getValue(q); })
          return static_cast<const Derived&>(*this).getValue(ip);
        else
          return static_cast<const Derived&>(*this).getValue(ip.getPoint());
      }

      /// @brief Returns an optional polynomial order bound.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(poly);
      }

      virtual ComplexFunctionBase* copy() const noexcept override = 0;
  };

  template <>
  class ComplexFunction<Integer> final
    : public ComplexFunctionBase<ComplexFunction<Integer>>
  {
    public:
      /// @brief Parent class type.
      using Parent = ComplexFunctionBase<ComplexFunction<Integer>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs a ComplexFunction from an integer value.
       * @param[in] x Constant integer value
       */
      ComplexFunction(const Integer& x)
        : m_x(x)
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      /// @brief Returns the stored constant value.
      constexpr
      const Integer& getValue() const
      {
        return m_x;
      }

      /// @brief Evaluates the constant function at a point.
      constexpr
      Complex getValue(const Geometry::Point&) const
      {
        return Complex(m_x, 0);
      }

      /// @brief Returns zero polynomial order for a constant function.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      const Integer m_x;
  };

  /// @brief Deduction guide for integer constants.
  ComplexFunction(Integer) -> ComplexFunction<Integer>;

  template <>
  class ComplexFunction<Real> final
    : public ComplexFunctionBase<ComplexFunction<Real>>
  {
    public:
      /// @brief Parent class type.
      using Parent = ComplexFunctionBase<ComplexFunction<Real>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs a ComplexFunction from an Real value.
       * @param[in] x Constant Real value
       */
      ComplexFunction(const Real& x)
        : m_x(x)
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      /// @brief Returns the stored constant value.
      constexpr
      const Real& getValue() const
      {
        return m_x;
      }

      /// @brief Evaluates the constant function at a point.
      constexpr
      Complex getValue(const Geometry::Point&) const
      {
        return Complex(m_x, 0);
      }

      /// @brief Returns zero polynomial order for a constant function.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      const Real m_x;
  };

  /// @brief Deduction guide for real constants.
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
      /// @brief Parent class type.
      using Parent = ComplexFunctionBase<ComplexFunction<Complex>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs a ComplexFunction from a constant scalar value.
       * @param[in] x Constant scalar value
       */
      ComplexFunction(const Complex& x)
        : m_x(x)
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_x(other.m_x)
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_x(other.m_x)
      {}

      /// @brief Returns the stored constant value.
      constexpr
      const Complex& getValue() const
      {
        return m_x;
      }

      /// @brief Evaluates the constant function at a point.
      constexpr
      Complex getValue(const Geometry::Point&) const
      {
        return m_x;
      }

      /// @brief Returns zero polynomial order for a constant function.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      ComplexFunction* copy() const noexcept override
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
   * @brief Complex function wrapper around a nested function.
   */
  template <class NestedDerived>
  class ComplexFunction<FunctionBase<NestedDerived>> final
    : public ComplexFunctionBase<ComplexFunction<NestedDerived>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Complex;

      /// @brief Parent class type.
      using Parent = ComplexFunctionBase<ComplexFunction<NestedDerived>>;

      using Parent::traceOf;

      using Parent::operator();

      /// @brief Constructs from a nested function.
      ComplexFunction(const FunctionBase<NestedDerived>& nested)
        : m_nested(nested.copy())
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_nested(other.m_nested->copy())
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_nested(std::move(other.m_nested))
      {}

      /// @brief Evaluates the nested function at a geometry point.
      constexpr
      ScalarType getValue(const Geometry::Point& v) const
      {
        return m_nested->getValue(v);
      }

      /// @brief Evaluates the nested function at an integration point.
      constexpr
      ScalarType getValue(const IntegrationPoint& ip) const
      {
        return m_nested->getValue(ip);
      }

      /// @brief Forwards trace-domain configuration to the nested function.
      template <class ... Args>
      constexpr
      ComplexFunction& traceOf(const Args&... args)
      {
        m_nested->traceOf(args...);
        return *this;
      }

      /// @brief Returns the nested function order.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return m_nested->getOrder(poly);
      }

      ComplexFunction* copy() const noexcept override
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
   * @brief Complex function built from real and imaginary functions.
   */
  template <class RealNestedDerived, class ImagNestedDerived>
  class ComplexFunction<FunctionBase<RealNestedDerived>, FunctionBase<ImagNestedDerived>> final
    : public ComplexFunctionBase<ComplexFunction<FunctionBase<RealNestedDerived>, FunctionBase<ImagNestedDerived>>>
  {
    public:
      /// @brief Real-part function type.
      using RealFunctionType =
        FunctionBase<RealNestedDerived>;

      /// @brief Imaginary-part function type.
      using ImagFunctionType =
        FunctionBase<ImagNestedDerived>;

      /// @brief Real-part function range type.
      using RealFunctionRangeType =
        typename FormLanguage::Traits<RealFunctionType>::RangeType;

      /// @brief Imaginary-part function range type.
      using ImagFunctionRangeType =
        typename FormLanguage::Traits<ImagFunctionType>::RangeType;

      /// @brief Parent class type.
      using Parent =
        ComplexFunctionBase<
          ComplexFunction<
            FunctionBase<RealNestedDerived>,
            FunctionBase<ImagNestedDerived>>>;

      using Parent::traceOf;

      using Parent::operator();

      static_assert(std::is_same_v<RealFunctionRangeType, Real>);

      static_assert(std::is_same_v<ImagFunctionRangeType, Real>);

      /// @brief Constructs from real and imaginary functions.
      ComplexFunction(const RealFunctionType& re, const ImagFunctionType& imag)
        : m_re(re.copy()), m_imag(imag.copy())
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_re(other.m_re->copy()), m_imag(other.m_imag->copy())
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_re(std::move(other.m_re)), m_imag(std::move(other.m_imag))
      {}

      /// @brief Evaluates the complex function at a geometry point.
      constexpr
      Complex getValue(const Geometry::Point& p) const
      {
        return { m_re->getValue(p), m_imag->getValue(p) };
      }

      /// @brief Evaluates the complex function at an integration point.
      constexpr
      Complex getValue(const IntegrationPoint& ip) const
      {
        Complex res;
        if constexpr (requires { m_re->getValue(ip); })
          res.real(m_re->getValue(ip));
        else
          res.real(m_re->getValue(ip.getPoint()));

        if constexpr (requires { m_imag->getValue(ip); })
          res.imag(m_imag->getValue(ip));
        else
          res.imag(m_imag->getValue(ip.getPoint()));
        return res;
      }

      /// @brief Forwards trace-domain configuration to both parts.
      template <class ... Args>
      constexpr
      ComplexFunction& traceOf(const Args&... args)
      {
        m_re->traceOf(args...);
        m_imag->traceOf(args...);
        return *this;
      }

      /// @brief Returns the maximum order of the real and imaginary parts.
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const auto reOrder = m_re->getOrder(geom);
        const auto imOrder = m_imag->getOrder(geom);

        if (!reOrder || !imOrder)
          return std::nullopt;

        return std::max(*reOrder, *imOrder);
      }

      ComplexFunction* copy() const noexcept override
      {
        return new ComplexFunction(*this);
      }

    private:
      std::unique_ptr<RealFunctionType> m_re;
      std::unique_ptr<ImagFunctionType> m_imag;
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
      /// @brief Parent class type.
      using Parent = ComplexFunctionBase<ComplexFunction<F>>;

      using Parent::traceOf;

      using Parent::operator();

      /// @brief Constructs from a callable.
      ComplexFunction(const F& f)
        : m_f(f)
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_f(other.m_f)
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f))
      {}

      /// @brief Evaluates the callable at a geometry point.
      constexpr
      Complex getValue(const Geometry::Point& v) const
      {
        return m_f(v);
      }

      /// @brief Evaluates the callable at an integration point.
      constexpr
      Complex getValue(const IntegrationPoint& ip) const
      {
        if constexpr (std::is_invocable_r_v<Complex, F, const IntegrationPoint&>)
          return m_f(ip);
        else
          return m_f(ip.getPoint());
      }

      /// @brief Returns no polynomial order for an arbitrary callable.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      ComplexFunction* copy() const noexcept override
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

  /// @brief Complex function built from real and imaginary callables.
  template <class FReal, class FImag>
  class ComplexFunction<FReal, FImag> final
    : public ComplexFunctionBase<ComplexFunction<FReal, FImag>>
  {
    static_assert(std::is_invocable_v<FReal, const Geometry::Point&>);
    static_assert(std::is_invocable_v<FImag, const Geometry::Point&>);

    public:
      /// @brief Parent class type.
      using Parent = ComplexFunctionBase<ComplexFunction<FReal, FImag>>;

      using Parent::traceOf;

      using Parent::operator();

      /// @brief Constructs from real and imaginary callables.
      ComplexFunction(const FReal& re, const FImag& imag)
        : m_re(re), m_imag(imag)
      {}

      /// @brief Copy constructor.
      ComplexFunction(const ComplexFunction& other)
        : Parent(other),
          m_re(other.m_re), m_imag(other.m_imag)
      {}

      /// @brief Move constructor.
      ComplexFunction(ComplexFunction&& other)
        : Parent(std::move(other)),
          m_re(std::move(other.m_re)), m_imag(std::move(other.m_imag))
      {}

      /// @brief Evaluates the complex function at a geometry point.
      constexpr
      Complex getValue(const Geometry::Point& p) const
      {
        return { m_re(p), m_imag(p) };
      }

      /// @brief Evaluates the complex function at an integration point.
      constexpr
      Complex getValue(const IntegrationPoint& ip) const
      {
        Complex res;
        if constexpr (std::is_invocable_v<FReal, const IntegrationPoint&>)
          res.real(m_re(ip));
        else
          res.real(m_re(ip.getPoint()));

        if constexpr (std::is_invocable_v<FImag, const IntegrationPoint&>)
          res.imag(m_imag(ip));
        else
          res.imag(m_imag(ip.getPoint()));
        return res;
      }

      /// @brief Returns no polynomial order for arbitrary callables.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      ComplexFunction* copy() const noexcept override
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

/// @endcond
#endif
