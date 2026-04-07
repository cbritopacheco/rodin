/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Common.h
 * @brief Common mathematical functions and operations.
 *
 * This file provides a comprehensive collection of mathematical functions including
 * elementary functions, trigonometric functions, combinatorics, and linear algebra
 * operations. Functions are templated to work with various numeric types.
 */
#ifndef RODIN_MATH_COMMON_H
#define RODIN_MATH_COMMON_H

#include <cmath>
#include <Eigen/Core>

#include "Rodin/Types.h"

#include "ForwardDecls.h"

namespace Rodin::Math
{
  /**
   * @brief Computes the absolute value of a value.
   *
   * Returns @f$ |x| @f$, the absolute value (magnitude) of @f$ x @f$.
   *
   * @tparam T Type of value (Real, Complex, etc.)
   * @param[in] x Value
   * @return Absolute value @f$ |x| @f$
   */
  template <class T>
  constexpr
  auto abs(const T& x)
  {
    return x < static_cast<T>(0.0) ? -x : x;
  }

  /**
   * @brief Computes the exponential function.
   *
   * Returns @f$ e^x @f$ where @f$ e \approx 2.71828 @f$ is Euler's number.
   *
   * @tparam T Type of value
   * @param[in] x Exponent
   * @return @f$ e^x @f$
   */
  template <class T>
  constexpr
  auto exp(const T& x)
  {
    return std::exp(x);
  }

  /**
   * @brief Computes the complex conjugate of a complex number.
   *
   * For @f$ z = a + bi @f$, returns @f$ \bar{z} = a - bi @f$.
   *
   * @param[in] x Complex number
   * @return Complex conjugate @f$ \bar{x} @f$
   */
  constexpr
  Complex conj(const Complex& x)
  {
    return std::conj(x);
  }

  /**
   * @brief Computes the conjugate of an Eigen matrix.
   *
   * For matrices with complex entries, returns the element-wise conjugate.
   *
   * @tparam T Eigen matrix type
   * @param[in] x Matrix
   * @return Conjugate matrix
   */
  template <class T>
  constexpr
  auto conj(const Eigen::MatrixBase<T>& x)
  {
    return x.conjugate();
  }

  /**
   * @brief Identity conjugate for real numbers.
   *
   * For real numbers, the conjugate is the number itself.
   *
   * @param[in] x Real number
   * @return @f$ x @f$
   */
  constexpr
  Real conj(const Real& x)
  {
    return x;
  }

  /**
   * @brief Computes the square of a value.
   *
   * Returns @f$ x^2 = x \times x @f$ efficiently without calling std::pow.
   *
   * @tparam Base Type of the value
   * @tparam Exponent Unused template parameter for compatibility
   * @param[in] base Value to square
   * @return @f$ \text{base}^2 @f$
   */
  template <class Base>
  constexpr
  auto pow2(const Base& base)
  {
    return base * base;
  }

  /**
   * @brief Computes the compile-time integer power.
   *
   * Returns @f$ x^N @f$ where @f$ N @f$ is known at compile time, using
   * exponentiation by squaring for efficiency.
   *
   * @tparam N Exponent (non-negative integer)
   * @tparam Base Type of the base value
   * @param[in] x Base value
   * @return @f$ x^N @f$
   */
  template <size_t N, class Base>
  constexpr Base pow(const Base& x)
  {
    if constexpr (N == 0)
    {
      return Base(1);
    }
    else if constexpr (N == 1)
    {
      return x;
    }
    else if constexpr (N % 2 == 0)
    {
      const Base y = pow<N / 2>(x);
      return y * y;
    }
    else
    {
      const Base y = pow<N / 2>(x);
      return y * y * x;
    }
  }

  /**
   * @brief Computes the compile-time integer power via integral_constant.
   *
   * Convenience overload accepting the exponent as a
   * `std::integral_constant<size_t, N>` parameter.
   *
   * @tparam Base Type of the base value
   * @tparam N Exponent
   * @param[in] x Base value
   * @return @f$ x^N @f$
   */
  template <class Base, size_t N>
  constexpr Base pow(const Base& x, std::integral_constant<size_t, N>)
  {
    return pow<N>(x);
  }

  /**
   * @brief Computes a power with arbitrary exponent.
   *
   * Returns @f$ x^y @f$ where @f$ x @f$ is the base and @f$ y @f$ is the exponent.
   *
   * @tparam Base Type of the base
   * @tparam Exponent Type of the exponent
   * @param[in] base Base value @f$ x @f$
   * @param[in] exponent Exponent value @f$ y @f$
   * @return @f$ x^y @f$
   */
  template <class Base, class Exponent>
  constexpr
  auto pow(const Base& base, const Exponent& exponent)
  {
    return std::pow(base, exponent);
  }

  /**
   * @brief Computes the square root of a value.
   *
   * Returns @f$ \sqrt{x} @f$, the principal square root of @f$ x @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value (must be non-negative for real types)
   * @return @f$ \sqrt{x} @f$
   */
  template <class T>
  constexpr
  auto sqrt(const T& x)
  {
    return std::sqrt(x);
  }

  /**
   * @brief Determines if a floating point number is not-a-number (NaN).
   *
   * Checks if @f$ x @f$ is NaN (not a number), which typically results from
   * undefined operations like @f$ 0/0 @f$ or @f$ \sqrt{-1} @f$ for reals.
   *
   * @tparam T Type of value
   * @param[in] x Value to check
   * @return True if @f$ x @f$ is NaN, false otherwise
   */
  template <class T>
  constexpr
  Boolean isNaN(const T& x)
  {
    return std::isnan(x);
  }

  /**
   * @brief Determines if a complex number is not-a-number (NaN).
   *
   * Returns true if either the real or imaginary part is NaN.
   *
   * @param[in] x Complex number to check
   * @return True if real or imaginary part is NaN, false otherwise
   */
  inline
  constexpr
  Boolean isNaN(const Complex& x)
  {
    return std::isnan(x.real()) || std::isnan(x.imag());
  }

  /**
   * @brief Determines if a floating point number is infinite.
   *
   * Checks if @f$ x = \pm\infty @f$, which can result from overflow or
   * operations like @f$ 1/0 @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value to check
   * @return True if @f$ x @f$ is @f$ +\infty @f$ or @f$ -\infty @f$, false otherwise
   */
  template <class T>
  constexpr
  Boolean isInf(const T& x)
  {
    return std::isinf(x);
  }

  /**
   * @brief Computes the cosine function.
   *
   * Returns @f$ \cos(x) @f$ where @f$ x @f$ is in radians.
   *
   * @tparam T Type of value
   * @param[in] x Angle in radians
   * @return @f$ \cos(x) @f$
   */
  template <class T>
  constexpr
  auto cos(const T& x)
  {
    return std::cos(x);
  }

  /**
   * @brief Computes the hyperbolic cosine function.
   *
   * Returns @f$ \cosh(x) = \frac{e^x + e^{-x}}{2} @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value
   * @return @f$ \cosh(x) @f$
   */
  template <class T>
  constexpr
  auto cosh(const T& x)
  {
    return std::cosh(x);
  }

  /**
   * @brief Computes the sine function.
   *
   * Returns @f$ \sin(x) @f$ where @f$ x @f$ is in radians.
   *
   * @tparam T Type of value
   * @param[in] x Angle in radians
   * @return @f$ \sin(x) @f$
   */
  template <class T>
  constexpr
  auto sin(const T& x)
  {
    return std::sin(x);
  }

  /**
   * @brief Computes the hyperbolic sine function.
   *
   * Returns @f$ \sinh(x) = \frac{e^x - e^{-x}}{2} @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value
   * @return @f$ \sinh(x) @f$
   */
  template <class T>
  constexpr
  auto sinh(const T& x)
  {
    return std::sinh(x);
  }

  /**
   * @brief Computes the tangent function.
   *
   * Returns @f$ \tan(x) = \frac{\sin(x)}{\cos(x)} @f$ where @f$ x @f$ is in radians.
   *
   * @tparam T Type of value
   * @param[in] x Angle in radians
   * @return @f$ \tan(x) @f$
   */
  template <class T>
  constexpr
  auto tan(const T& x)
  {
   return std::tan(x);
  }

  /**
   * @brief Computes the sign (signum) function.
   *
   * Returns:
   * @f[
   *   \text{sgn}(x) = \begin{cases}
   *     -1 & \text{if } x < 0 \\
   *     0 & \text{if } x = 0 \\
   *     +1 & \text{if } x > 0
   *   \end{cases}
   * @f]
   *
   * @tparam T Type of value
   * @param[in] x Value
   * @return Sign of @f$ x @f$: -1, 0, or +1
   */
  template <typename T>
  constexpr
  T sgn(const T& x)
  {
    return (T(0) < x) - (x < T(0));
  }

  /**
   * @brief Computes the binomial coefficient.
   *
   * Returns @f$ \binom{n}{k} = \frac{n!}{k!(n-k)!} @f$, the number of ways to
   * choose @f$ k @f$ items from @f$ n @f$ items without regard to order.
   *
   * @tparam T Type of value (typically Integer)
   * @param[in] n Total number of items
   * @param[in] k Number of items to choose
   * @return @f$ \binom{n}{k} @f$
   * @pre @f$ 0 \leq k \leq n @f$
   */
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

  /**
   * @brief Computes the factorial function.
   *
   * Returns @f$ n! = 1 \times 2 \times 3 \times \cdots \times n @f$.
   *
   * @tparam T Type of value (typically Integer)
   * @param[in] n Non-negative integer
   * @return @f$ n! @f$
   * @pre @f$ n \geq 0 @f$
   */
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

  /**
   * @brief Computes the permutation function.
   *
   * Returns @f$ P(n, k) = \frac{n!}{(n-k)!} @f$, the number of ways to arrange
   * @f$ k @f$ items from @f$ n @f$ items where order matters.
   *
   * @tparam T Type of value (typically Integer)
   * @param[in] n Total number of items
   * @param[in] k Number of items to arrange
   * @return @f$ P(n, k) @f$
   * @pre @f$ 0 \leq k \leq n @f$
   */
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

  /**
   * @brief Returns a quiet NaN (not-a-number) value.
   *
   * Returns a NaN value of the specified type, useful for indicating
   * undefined or invalid results.
   *
   * @tparam T Type of value (must be a floating point type)
   * @return Quiet NaN value of type T
   */
  template <class T>
  constexpr
  auto nan()
  {
    if constexpr (std::is_same_v<T, Complex>)
    {
      return Complex(std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN());
    }
    else
    {
      return std::numeric_limits<T>::quiet_NaN();
    }
  }

  /**
   * @brief Computes the sum of two values.
   *
   * Returns @f$ \text{lhs} + \text{rhs} @f$. This function is used internally
   * for type deduction in the form language.
   *
   * @tparam LHS Type of left-hand side
   * @tparam RHS Type of right-hand side
   * @param[in] lhs Left-hand side operand
   * @param[in] rhs Right-hand side operand
   * @return @f$ \text{lhs} + \text{rhs} @f$
   */
  template <class LHS, class RHS>
  constexpr
  auto sum(const LHS& lhs, const RHS& rhs)
  {
    return lhs + rhs;
  }

  /**
   * @brief Computes the unary minus (negation).
   *
   * Returns @f$ -\text{op} @f$. This function is used internally for type
   * deduction in the form language.
   *
   * @tparam Operand Type of operand
   * @param[in] op Operand to negate
   * @return @f$ -\text{op} @f$
   */
  template <class Operand>
  constexpr
  auto minus(const Operand& op)
  {
    return -op;
  }

  /**
   * @brief Computes the difference of two values.
   *
   * Returns @f$ \text{lhs} - \text{rhs} @f$. This function is used internally
   * for type deduction in the form language.
   *
   * @tparam LHS Type of left-hand side
   * @tparam RHS Type of right-hand side
   * @param[in] lhs Left-hand side operand
   * @param[in] rhs Right-hand side operand
   * @return @f$ \text{lhs} - \text{rhs} @f$
   */
  template <class LHS, class RHS>
  constexpr
  auto minus(const LHS& lhs, const RHS& rhs)
  {
    return lhs - rhs;
  }

  /**
   * @brief Computes the product of two values.
   *
   * Returns @f$ \text{lhs} \times \text{rhs} @f$. This function is used internally
   * for type deduction in the form language.
   *
   * @tparam LHS Type of left-hand side
   * @tparam RHS Type of right-hand side
   * @param[in] lhs Left-hand side operand
   * @param[in] rhs Right-hand side operand
   * @return @f$ \text{lhs} \times \text{rhs} @f$
   */
  template <class LHS, class RHS>
  constexpr
  auto mult(const LHS& lhs, const RHS& rhs)
  {
    return lhs * rhs;
  }

  /**
   * @brief Computes the division of two values.
   *
   * Returns @f$ \text{lhs} / \text{rhs} @f$. This function is used internally
   * for type deduction in the form language.
   *
   * @tparam LHS Type of left-hand side
   * @tparam RHS Type of right-hand side
   * @param[in] lhs Left-hand side operand (numerator)
   * @param[in] rhs Right-hand side operand (denominator)
   * @return @f$ \text{lhs} / \text{rhs} @f$
   */
  template <class LHS, class RHS>
  constexpr
  auto division(const LHS& lhs, const RHS& rhs)
  {
    return lhs / rhs;
  }

  /**
   * @brief Computes the dot product of two real numbers.
   *
   * For real numbers, the dot product is simply multiplication:
   * @f$ \text{lhs} \cdot \text{rhs} = \text{lhs} \times \text{rhs} @f$.
   *
   * @param[in] lhs Left-hand side operand
   * @param[in] rhs Right-hand side operand
   * @return @f$ \text{lhs} \times \text{rhs} @f$
   */
  inline
  constexpr
  Real dot(const Real& lhs, const Real& rhs)
  {
    return lhs * rhs;
  }

  /**
   * @brief Computes the dot product of two complex numbers.
   *
   * For complex numbers, the dot product uses conjugation:
   * @f$ z_1 \cdot z_2 = z_1 \times \overline{z_2} @f$.
   *
   * @param[in] lhs Left-hand side complex number
   * @param[in] rhs Right-hand side complex number
   * @return @f$ \text{lhs} \times \overline{\text{rhs}} @f$
   */
  inline
  constexpr
  Complex dot(const Complex& lhs, const Complex& rhs)
  {
    return lhs * conj(rhs);
  }

  /**
   * @brief Computes the dot product of two Eigen vectors/matrices.
   *
   * For vectors @f$ \mathbf{u} @f$ and @f$ \mathbf{v} @f$, computes:
   * @f[
   *   \mathbf{u} \cdot \mathbf{v} = \sum_i u_i \overline{v_i}
   * @f]
   * where @f$ \overline{v_i} @f$ is the complex conjugate (identity for real values).
   *
   * @tparam LHSDerived Eigen type of left-hand side
   * @tparam RHSDerived Eigen type of right-hand side
   * @param[in] lhs Left-hand side vector/matrix
   * @param[in] rhs Right-hand side vector/matrix
   * @return Dot product (scalar)
   */
  template <class LHSDerived, class RHSDerived>
  constexpr
  auto dot(const Eigen::MatrixBase<LHSDerived>& lhs, const Eigen::MatrixBase<RHSDerived>& rhs)
  {
    assert(lhs.rows() == rhs.rows());
    assert(lhs.cols() == rhs.cols());
    return (lhs.array() * rhs.conjugate().array()).sum();
  }

  /**
   * @brief Computes the dot product of two SpatialVector objects.
   *
   * @tparam LHSScalar Scalar type of left-hand side
   * @tparam RHSScalar Scalar type of right-hand side
   * @param[in] lhs Left-hand side SpatialVector
   * @param[in] rhs Right-hand side SpatialVector
   * @return Dot product (scalar)
   */
  template <class LHSScalar, class RHSScalar>
  constexpr
  auto dot(const SpatialVector<LHSScalar>& lhs, const SpatialVector<RHSScalar>& rhs)
  {
    return lhs.dot(rhs);
  }

  /**
   * @brief Computes the dot product of two SpatialMatrix objects.
   *
   * @tparam LHSScalar Scalar type of left-hand side
   * @tparam RHSScalar Scalar type of right-hand side
   * @param[in] lhs Left-hand side SpatialMatrix
   * @param[in] rhs Right-hand side SpatialMatrix
   * @return Dot product (scalar)
   */
  template <class LHSScalar, class RHSScalar>
  constexpr
  auto dot(const SpatialMatrix<LHSScalar>& lhs, const SpatialMatrix<RHSScalar>& rhs)
  {
    return lhs.dot(rhs);
  }

  /**
   * @brief Computes the dot product of an Eigen vector and a SpatialVector.
   *
   * @tparam LHSDerived Eigen type of left-hand side
   * @tparam RHSScalar Scalar type of SpatialVector
   * @param[in] lhs Left-hand side Eigen vector
   * @param[in] rhs Right-hand side SpatialVector
   * @return Dot product (scalar)
   */
  template <class LHSDerived, class RHSScalar>
  constexpr
  auto dot(const Eigen::MatrixBase<LHSDerived>& lhs, const Math::SpatialVector<RHSScalar>& rhs)
  {
    return lhs.dot(rhs.getData().head(static_cast<Eigen::Index>(rhs.size())));
  }

  /**
   * @brief Computes the dot product of a SpatialVector and an Eigen vector.
   *
   * @tparam LHSScalar Scalar type of SpatialVector
   * @tparam RHSDerived Eigen type of right-hand side
   * @param[in] lhs Left-hand side SpatialVector
   * @param[in] rhs Right-hand side Eigen vector
   * @return Dot product (scalar)
   */
  template <class LHSScalar, class RHSDerived>
  constexpr
  auto dot(const Math::SpatialVector<LHSScalar>& lhs, const Eigen::MatrixBase<RHSDerived>& rhs)
  {
    return lhs.getData().head(static_cast<Eigen::Index>(lhs.size())).dot(rhs);
  }

  /**
   * @brief Computes the dot product of an Eigen matrix and a SpatialMatrix.
   *
   * @tparam LHSDerived Eigen type of left-hand side
   * @tparam RHSScalar Scalar type of SpatialMatrix
   * @param[in] lhs Left-hand side Eigen matrix
   * @param[in] rhs Right-hand side SpatialMatrix
   * @return Dot product (scalar)
   */
  template <class LHSDerived, class RHSScalar>
  constexpr
  auto dot(const Eigen::MatrixBase<LHSDerived>& lhs, const Math::SpatialMatrix<RHSScalar>& rhs)
  {
    return lhs.dot(rhs.getData().head(static_cast<Eigen::Index>(rhs.size())));
  }

  /**
   * @brief Computes the dot product of a SpatialMatrix and an Eigen matrix.
   *
   * @tparam LHSScalar Scalar type of SpatialMatrix
   * @tparam RHSDerived Eigen type of right-hand side
   * @param[in] lhs Left-hand side SpatialMatrix
   * @param[in] rhs Right-hand side Eigen matrix
   * @return Dot product (scalar)
   */
  template <class LHSScalar, class RHSDerived>
  constexpr
  auto dot(const Math::SpatialMatrix<LHSScalar>& lhs, const Eigen::MatrixBase<RHSDerived>& rhs)
  {
    return lhs.getData().head(static_cast<Eigen::Index>(lhs.size())).dot(rhs);
  }

  /**
   * @brief Computes the natural logarithm.
   *
   * Returns @f$ \ln(x) @f$, the natural (base-@f$ e @f$) logarithm of @f$ x @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value (must be positive for real types)
   * @return @f$ \ln(x) @f$
   */
  template <class T>
  constexpr
  auto log(const T& x)
  {
    return std::log(x);
  }

  /**
   * @brief Computes the base-2 logarithm.
   *
   * Returns @f$ \log_2(x) @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value (must be positive for real types)
   * @return @f$ \log_2(x) @f$
   */
  template <class T>
  constexpr
  auto log2(const T& x)
  {
    return std::log2(x);
  }

  /**
   * @brief Computes the base-10 logarithm.
   *
   * Returns @f$ \log_{10}(x) @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value (must be positive for real types)
   * @return @f$ \log_{10}(x) @f$
   */
  template <class T>
  constexpr
  auto log10(const T& x)
  {
    return std::log10(x);
  }

  /**
   * @brief Computes the arc cosine (inverse cosine).
   *
   * Returns @f$ \arccos(x) @f$ in radians, in the range @f$ [0, \pi] @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value in @f$ [-1, 1] @f$
   * @return @f$ \arccos(x) @f$ in radians
   */
  template <class T>
  constexpr
  auto acos(const T& x)
  {
    return std::acos(x);
  }

  /**
   * @brief Computes the arc sine (inverse sine).
   *
   * Returns @f$ \arcsin(x) @f$ in radians, in the range @f$ [-\pi/2, \pi/2] @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value in @f$ [-1, 1] @f$
   * @return @f$ \arcsin(x) @f$ in radians
   */
  template <class T>
  constexpr
  auto asin(const T& x)
  {
    return std::asin(x);
  }

  /**
   * @brief Computes the arc tangent (inverse tangent).
   *
   * Returns @f$ \arctan(x) @f$ in radians, in the range @f$ (-\pi/2, \pi/2) @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value
   * @return @f$ \arctan(x) @f$ in radians
   */
  template <class T>
  constexpr
  auto atan(const T& x)
  {
    return std::atan(x);
  }

  /**
   * @brief Computes the two-argument arc tangent.
   *
   * Returns @f$ \text{atan2}(y, x) @f$ in radians, in the range @f$ (-\pi, \pi] @f$.
   * This is the angle of the vector @f$ (x, y) @f$ from the positive @f$ x @f$-axis.
   *
   * @tparam T Type of value
   * @param[in] y The y-coordinate
   * @param[in] x The x-coordinate
   * @return @f$ \text{atan2}(y, x) @f$ in radians
   */
  template <class T>
  constexpr
  auto atan2(const T& y, const T& x)
  {
    return std::atan2(y, x);
  }

  /**
   * @brief Computes the hyperbolic tangent function.
   *
   * Returns @f$ \tanh(x) = \frac{\sinh(x)}{\cosh(x)} @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value
   * @return @f$ \tanh(x) @f$
   */
  template <class T>
  constexpr
  auto tanh(const T& x)
  {
    return std::tanh(x);
  }

  /**
   * @brief Returns the minimum of two values.
   *
   * @tparam T Type of value
   * @param[in] a First value
   * @param[in] b Second value
   * @return The smaller of @f$ a @f$ and @f$ b @f$
   */
  template <class T>
  constexpr
  auto min(const T& a, const T& b)
  {
    return (a < b) ? a : b;
  }

  /**
   * @brief Returns the maximum of two values.
   *
   * @tparam T Type of value
   * @param[in] a First value
   * @param[in] b Second value
   * @return The larger of @f$ a @f$ and @f$ b @f$
   */
  template <class T>
  constexpr
  auto max(const T& a, const T& b)
  {
    return (a > b) ? a : b;
  }

  /**
   * @brief Clamps a value to the range @f$ [\text{lo}, \text{hi}] @f$.
   *
   * If @f$ x < \text{lo} @f$, returns @f$ \text{lo} @f$.
   * If @f$ x > \text{hi} @f$, returns @f$ \text{hi} @f$.
   * Otherwise returns @f$ x @f$.
   *
   * @tparam T Type of value
   * @param[in] x Value to clamp
   * @param[in] lo Lower bound
   * @param[in] hi Upper bound
   * @return Clamped value in @f$ [\text{lo}, \text{hi}] @f$
   */
  template <class T>
  constexpr
  auto clamp(const T& x, const T& lo, const T& hi)
  {
    return min(max(x, lo), hi);
  }
}

#endif
