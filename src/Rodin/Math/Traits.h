/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_TRAITS_H
#define RODIN_MATH_TRAITS_H

#include <type_traits>

#include <Eigen/Core>

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Types.h"

#include "Common.h"

namespace Rodin::FormLanguage
{
  template <class T>
  struct IsEigenObject
  {
    static constexpr bool Value =
      std::is_base_of_v<Eigen::EigenBase<typename std::decay<T>::type>, typename std::decay<T>::type>;
  };

  template <>
  struct Traits<Boolean>
  {
    using ScalarType = Boolean;
  };

  template <>
  struct Traits<Integer>
  {
    using ScalarType = Integer;
  };

  template <>
  struct Traits<Real>
  {
    using ScalarType = Real;
  };

  template <>
  struct Traits<Complex>
  {
    using ScalarType = Complex;
  };

  template <class LHS, class RHS>
  struct Sum
  {
    using Type = decltype(Math::sum(std::declval<LHS>(), std::declval<RHS>()));
  };

  template <class LHS, class RHS>
  struct Minus
  {
    using Type = decltype(Math::minus(std::declval<LHS>(), std::declval<RHS>()));
  };

  template <class Operand>
  struct UnaryMinus
  {
    using Type = decltype(Math::minus(std::declval<Operand>()));
  };

  template <class LHS, class RHS>
  struct Mult
  {
    using Type = decltype(Math::mult(std::declval<LHS>(), std::declval<RHS>()));
  };

  template <class LHS, class RHS>
  struct Division
  {
    using Type = decltype(Math::division(std::declval<LHS>(), std::declval<RHS>()));
  };

  template <class LHS, class RHS>
  struct Dot
  {
    using Type = decltype(Math::dot(std::declval<LHS>(), std::declval<RHS>()));
  };
}

#endif

