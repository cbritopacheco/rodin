/**
 * @file Minus.h
 * @brief Subtraction operations for functions, shape functions, and integrators.
 *
 * Provides operator- overloads for computing differences between functions,
 * shape functions, and integrators. All subtraction operations are implemented
 * as Sum(lhs, UnaryMinus(rhs)) to reuse addition and negation logic.
 */

#ifndef RODIN_VARIATIONAL_MINUS_H
#define RODIN_VARIATIONAL_MINUS_H

#include <type_traits>

#include "ForwardDecls.h"
#include "Sum.h"
#include "UnaryMinus.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
  /**
   * @brief Subtracts a number from a function.
   * @f[
   *   (f - c)(x) = f(x) - c
   * @f]
   * @param lhs Function operand
   * @param rhs Number operand
   * @returns Difference function @f$ f - c @f$
   */
  template <class LHSDerived, class Number, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator-(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Sum(lhs, UnaryMinus(RealFunction(rhs)));
  }

  /**
   * @brief Subtracts a function from a number.
   * @f[
   *   (c - f)(x) = c - f(x)
   * @f]
   * @param lhs Number operand
   * @param rhs Function operand
   * @returns Difference function @f$ c - f @f$
   */
  template <class Number, class RHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator-(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(RealFunction(lhs), UnaryMinus(rhs));
  }

  /**
   * @brief Subtracts two functions.
   * @f[
   *   (f - g)(x) = f(x) - g(x)
   * @f]
   * @param lhs First function
   * @param rhs Second function
   * @returns Difference function @f$ f - g @f$
   */
  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator-(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }

  /**
   * @brief Subtracts two shape functions.
   * @f[
   *   (u - v)_i = u_i - v_i
   * @f]
   * @param lhs First shape function
   * @param rhs Second shape function
   * @returns Difference shape function @f$ u - v @f$
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  constexpr
  auto
  operator-(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs, const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }

  /**
   * @brief Subtracts two bilinear form integrators.
   * @param lhs First integrator
   * @param rhs Second integrator
   * @returns Difference integrator
   */
  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator-(
      const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
      const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }

  /**
   * @brief Subtracts a list of bilinear form integrators from a single integrator.
   * @param lhs Single integrator
   * @param rhs List of integrators
   * @returns Difference integrator list
   */
  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator-(
      const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }

  /**
   * @brief Subtracts a single bilinear form integrator from a list.
   * @param lhs List of integrators
   * @param rhs Single integrator
   * @returns Difference integrator list
   */
  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lhs,
      const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }

  /**
   * @brief Subtracts two lists of bilinear form integrators.
   * @param lhs First integrator list
   * @param rhs Second integrator list
   * @returns Difference integrator list
   */
  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }
}

#endif
