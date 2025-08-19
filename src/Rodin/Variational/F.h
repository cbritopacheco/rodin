/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_F_H
#define RODIN_VARIATIONAL_F_H

#include "RealFunction.h"

/**
 * @brief Contains built-in functions.
 */
namespace Rodin::Variational::F
{
  class X : public RealFunctionBase<X>
  {
    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<X>;

      X() = default;

      X(const X& other)
        : Parent(other)
      {}

      X(X&& other)
        : Parent(std::move(other))
      {}

      ScalarType getValue(const Geometry::Point& p) const
      {
        return p.x();
      }

      X* copy() const noexcept override
      {
        return new X(*this);
      }
  };

  /**
   * @brief Represents the first coordinate of the point.
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \cdots, x_d) = x_1 \: .
   * @f]
   */
  static const X x;

  class Y : public RealFunctionBase<Y>
  {
    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<Y>;

      Y() = default;

      Y(const Y& other)
        : Parent(other)
      {}

      Y(Y&& other)
        : Parent(std::move(other))
      {}

      ScalarType getValue(const Geometry::Point& p) const
      {
        return p.y();
      }

      Y* copy() const noexcept override
      {
        return new Y(*this);
      }
  };

  /**
   * @brief Represents the first coordinate of the point.
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \cdots, x_d) = x_2 \: .
   * @f]
   */
  static const Y y;

  class Z : public RealFunctionBase<Z>
  {
    public:
      using ScalarType = Real;

      using Parent = RealFunctionBase<Z>;

      Z() = default;

      Z(const Z& other)
        : Parent(other)
      {}

      Z(Z&& other)
        : Parent(std::move(other))
      {}

      ScalarType getValue(const Geometry::Point& p) const
      {
        return p.z();
      }

      Z* copy() const noexcept override
      {
        return new Z(*this);
      }
  };

  /**
   * @brief Represents the first coordinate of the point.
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \cdots, x_d) = x_3 \: .
   * @f]
   */
  static const Z z;
}

#endif
