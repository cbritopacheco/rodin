/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_F_H
#define RODIN_VARIATIONAL_F_H

/**
 * @file
 * @brief Coordinate projection functions for variational formulations
 *
 * This file defines built-in coordinate projection functions (x, y, z) that
 * extract individual components from geometric points. These are commonly used
 * in defining boundary conditions, source terms, and exact solutions.
 */

#include "RealFunction.h"

/**
 * @brief Contains built-in coordinate functions
 *
 * This namespace provides predefined real-valued functions that project
 * spatial coordinates, enabling convenient construction of position-dependent
 * expressions in variational formulations.
 */
namespace Rodin::Variational::F
{
  /**
   * @brief First coordinate projection function
   *
   * Class representing the @f$ x_1 @f$ coordinate projection.
   */
  class X : public RealFunctionBase<X>
  {
    public:
      /// @brief Parent class type
      using Parent = RealFunctionBase<X>;

      /**
       * @brief Default constructor
       */
      X() = default;

      /**
       * @brief Copy constructor
       * @param other X function to copy from
       */
      X(const X& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       * @param other X function to move from
       */
      X(X&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Evaluates the x-coordinate at a point
       * @param p Point at which to evaluate
       * @returns The x-coordinate of the point
       */
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return p.x();
      }

      /**
       * @brief Creates a copy of this function
       * @returns Pointer to newly allocated copy
       */
      X* copy() const noexcept override
      {
        return new X(*this);
      }
  };

  /**
   * @brief First coordinate projection @f$ x_1 @f$
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \ldots, x_d) = x_1
   * @f]
   *
   * @par Usage Example:
   * @code{.cpp}
   * using namespace Rodin::Variational;
   * auto bc = DirichletBC(u, F::x); // u = x on boundary
   * @endcode
   */
  static const X x;

  /**
   * @brief Second coordinate projection function
   *
   * Class representing the @f$ x_2 @f$ coordinate projection.
   */
  class Y : public RealFunctionBase<Y>
  {
    public:
      /// @brief Parent class type
      using Parent = RealFunctionBase<Y>;

      /**
       * @brief Default constructor
       */
      Y() = default;

      /**
       * @brief Copy constructor
       * @param other Y function to copy from
       */
      Y(const Y& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       * @param other Y function to move from
       */
      Y(Y&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Evaluates the y-coordinate at a point
       * @param p Point at which to evaluate
       * @returns The y-coordinate of the point
       */
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return p.y();
      }

      /**
       * @brief Creates a copy of this function
       * @returns Pointer to newly allocated copy
       */
      Y* copy() const noexcept override
      {
        return new Y(*this);
      }
  };

  /**
   * @brief Second coordinate projection @f$ x_2 @f$
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \ldots, x_d) = x_2
   * @f]
   *
   * @par Usage Example:
   * @code{.cpp}
   * using namespace Rodin::Variational;
   * auto bc = DirichletBC(u, F::y); // u = y on boundary
   * @endcode
   */
  static const Y y;

  /**
   * @brief Third coordinate projection function
   *
   * Class representing the @f$ x_3 @f$ coordinate projection.
   */
  class Z : public RealFunctionBase<Z>
  {
    public:
      /// @brief Parent class type
      using Parent = RealFunctionBase<Z>;

      /**
       * @brief Default constructor
       */
      Z() = default;

      /**
       * @brief Copy constructor
       * @param other Z function to copy from
       */
      Z(const Z& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       * @param other Z function to move from
       */
      Z(Z&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Evaluates the z-coordinate at a point
       * @param p Point at which to evaluate
       * @returns The z-coordinate of the point
       */
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return p.z();
      }

      /**
       * @brief Creates a copy of this function
       * @returns Pointer to newly allocated copy
       */
      Z* copy() const noexcept override
      {
        return new Z(*this);
      }
  };

  /**
   * @brief Third coordinate projection @f$ x_3 @f$
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \ldots, x_d) = x_3
   * @f]
   *
   * @par Usage Example:
   * @code{.cpp}
   * using namespace Rodin::Variational;
   * auto bc = DirichletBC(u, F::z); // u = z on boundary
   * @endcode
   */
  static const Z z;
}

#endif
