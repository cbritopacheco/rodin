/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_COMMON_H
#define RODIN_CORE_COMMON_H

#include <cmath>
#include <type_traits>
#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::Math
{
  /**
   * @brief Computes the absolute value of a value of type T.
   * @param[in] x Value
   * @tparam T Type of value
   * @returns Absolute of value
   */
  template <class T>
  constexpr
  auto abs(const T& x)
  {
    return std::abs(x);
  }

  template <class T>
  constexpr
  auto exp(const T& x)
  {
    return std::exp(x);
  }

  constexpr
  Complex conj(const Complex& x)
  {
    return std::conj(x);
  }

  template <class T>
  constexpr
  auto conj(const Eigen::MatrixBase<T>& x)
  {
    return x.conjugate();
  }

  constexpr
  Real conj(const Real& x)
  {
    return x;
  }

  template <class Base, class Exponent>
  constexpr
  auto pow2(const Base& base)
  {
    return base * base;
  }

  template <class Base, class Exponent>
  constexpr
  auto pow(const Base& base, const Exponent& exponent)
  {
    return std::pow(base, exponent);
  }

  /**
  * @brief Computes the square root of a value of type T.
  * @param[in] x Value
  * @tparam T Type of value
  * @returns Square root of value
  */
  template <class T>
  constexpr
  auto sqrt(const T& x)
  {
    return std::sqrt(x);
  }

  /**
  * @brief Determines if the floating point number is not-a-number (NaN).
  * @param[in] x Value
  * @tparam T Type of value
  * @returns True if value is NaN, false otherwise.
  */
  template <class T>
  constexpr
  Boolean isNaN(const T& x)
  {
    return std::isnan(x);
  }

  /**
  * @brief Determines if the floating point number is positive or negative infinity.
  * @param[in] x Value
  * @tparam T Type of value
  * @returns True if value is Inf (or -Inf), false otherwise.
  */
  template <class T>
  constexpr
  Boolean isInf(const T& x)
  {
    return std::isinf(x);
  }

  template <class T>
  constexpr
  auto cos(const T& x)
  {
    return std::cos(x);
  }

  template <class T>
  constexpr
  auto cosh(const T& x)
  {
    return std::cosh(x);
  }

  template <class T>
  constexpr
  auto sin(const T& x)
  {
    return std::sin(x);
  }

  template <class T>
  constexpr
  auto sinh(const T& x)
  {
    return std::sinh(x);
  }

  template <class T>
  constexpr
  auto tan(const T& x)
  {
   return std::tan(x);
  }

  template <typename T>
  constexpr
  T sgn(const T& x)
  {
    return (T(0) < x) - (x < T(0));
  }

  template <class T>
  constexpr
  T binom(const T& n, const T& k)
  {
    assert(T(0) <= n);
    assert(T(0) <= k);
    assert(k <= n);
    T res(1);
    for (T i = 0; i < std::min(k, n - k); ++i)
    {
      res *= (n - i);
      res /= (i + T(1));
    }
    return res;
  }

  template <class T>
  constexpr
  T factorial(const T& n)
  {
    assert(T(0) <= n);
    T res(1);
    for (T i = T(2); i <= n; ++i)
      res *= i;
    return res;
  }

  template <class T>
  constexpr
  T permutation(const T& n, const T& k)
  {
    assert(T(0) <= n);
    assert(T(0) <= k);
    assert(k <= n);
    T res(1);
    for (T i = 0; i < k; i++)
        res *= (n - i);
    return res;
  }

  template <class T>
  constexpr
  auto nan()
  {
    return std::numeric_limits<T>::quiet_NaN();
  }

  template <class LHS, class RHS>
  constexpr
  auto sum(const LHS& lhs, const RHS& rhs)
  {
    return lhs + rhs;
  }

  template <class Operand>
  constexpr
  auto minus(const Operand& op)
  {
    return -op;
  }

  template <class LHS, class RHS>
  constexpr
  auto minus(const LHS& lhs, const RHS& rhs)
  {
    return lhs - rhs;
  }

  template <class LHS, class RHS>
  constexpr
  auto mult(const LHS& lhs, const RHS& rhs)
  {
    return lhs * rhs;
  }

  template <class LHS, class RHS>
  constexpr
  auto division(const LHS& lhs, const RHS& rhs)
  {
    return lhs / rhs;
  }

  constexpr
  Real dot(const Real& lhs, const Real& rhs)
  {
    return lhs * rhs;
  }

  constexpr
  Complex dot(const Complex& lhs, const Complex& rhs)
  {
    return conj(lhs) * rhs;
  }

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto dot(const Eigen::MatrixBase<LHSDerived>& lhs, const Eigen::MatrixBase<RHSDerived>& rhs)
  {
    using LHS = Eigen::MatrixBase<LHSDerived>;
    using RHS = Eigen::MatrixBase<RHSDerived>;
    if constexpr (LHS::IsVectorAtCompileTime)
    {
      static_assert(RHS::IsVectorAtCompileTime);
      return lhs.dot(rhs);
    }
    else
    {
      return (lhs.conjugate().array() * rhs.array()).rowwise().sum().colwise().sum().value();
    }
  }
}

namespace Rodin::FormLanguage
{
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
