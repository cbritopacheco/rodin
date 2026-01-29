/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SIN_H
#define RODIN_VARIATIONAL_SIN_H

/**
 * @file Sine.h
 * @brief Sine trigonometric function operator for scalar functions.
 */

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup SinSpecializations Sin Template Specializations
   * @brief Template specializations of the Sin class.
   * @see Sin
   */

  /**
   * @brief Sine function operator for real-valued scalar functions.
   *
   * Applies the sine function pointwise to a given function:
   * @f[
   *    \text{Sin}(f)(x) = \sin(f(x))
   * @f]
   *
   * Common applications include:
   * - Periodic boundary conditions
   * - Wave equations and vibration problems
   * - Fourier series representations
   * - Trigonometric manufactured solutions
   *
   * @note Always returns a real-valued function for real input.
   * @see Cos, Tan, Sinh
   */

  /**
   * @ingroup SinSpecializations
   */
  template <class NestedDerived>
  class Sin<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Sin<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Sin<FunctionBase<NestedDerived>>>;

      /**
       * @brief Constructs sine of a function.
       * @param v Function to apply sine to
       */
      Sin(const OperandType& v)
        : m_operand(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Sin object to copy from
       */
      Sin(const Sin& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Sin object to move from
       */
      Sin(Sin&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Restricts operand to a single attribute.
       * @param attr Mesh attribute for trace restriction
       * @returns Reference to this object
       */
      constexpr
      Sin& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      /**
       * @brief Restricts operand to multiple attributes.
       * @param attrs Set of mesh attributes for trace restriction
       * @returns Reference to this object
       */
      constexpr
      Sin& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      /**
       * @brief Evaluates sine at a point.
       * @param p Point at which to evaluate
       * @returns @f$ \sin(f(p)) @f$
       */
      Real getValue(const Geometry::Point& p) const
      {
        return Math::sin(getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      Optional<size_t> getOrder(const Geometry::Polytope& g) const
      {
        const auto o = getOperand().getOrder(g);

        // Only constant -> constant preserves polynomial nature.
        if (o && *o == 0)
          return size_t{0};

        // Anything else: not a polynomial integrand in general.
        return {};
      }

      Sin* copy() const noexcept override
      {
        return new Sin(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Sin(const FunctionBase<NestedDerived>&) -> Sin<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Sin.
   */
  template <class NestedDerived>
  auto sin(const FunctionBase<NestedDerived>& f)
  {
    return Sin(f);
  }
}

#endif

