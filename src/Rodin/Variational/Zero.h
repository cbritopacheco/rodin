/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ZERO_H
#define RODIN_VARIATIONAL_ZERO_H

/**
 * @file
 * @brief Zero function for variational formulations
 *
 * This file defines the Zero class template, which represents the constant
 * zero function in variational problems. It provides both scalar and vector
 * specializations for use in initial conditions, boundary values, and source terms.
 */

#include <Rodin/Math/Common.h>

#include "VectorFunction.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ZeroSpecializations Zero Template Specializations
   * @brief Template specializations of the Zero class.
   * @see Zero
   */

  /**
   * @brief Scalar zero function
   * @tparam Scalar Scalar type (e.g., Real, Complex)
   * @ingroup ZeroSpecializations
   *
   * Represents the constant zero scalar function:
   * @f[
   *   f(x) = 0 \quad \forall x \in \Omega
   * @f]
   *
   * @see Zero
   */
  template <class Scalar>
  class Zero<Scalar> final
    : public ScalarFunctionBase<Scalar, Zero<Scalar>>
  {
    public:
      /// @brief Scalar value type
      using ScalarType = Scalar;

      /// @brief Parent class type
      using Parent = ScalarFunctionBase<Scalar, Zero<Scalar>>;

      /**
       * @brief Default constructor
       */
      Zero() {}

      /**
       * @brief Copy constructor
       * @param other Zero function to copy from
       */
      Zero(const Zero& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       * @param other Zero function to move from
       */
      Zero(Zero&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Evaluates the zero function at a point
       * @param[in] p Point at which to evaluate (unused)
       * @returns Always returns 0
       */
      constexpr
      ScalarType getValue(const Geometry::Point&) const
      {
        return 0;
      }

      /**
       * @brief Creates a copy of this zero function
       * @returns Pointer to newly allocated copy
       */
      Zero* copy() const noexcept override
      {
        return new Zero(*this);
      }
  };

  /**
   * @brief CTAD (Class Template Argument Deduction) for scalar Zero
   *
   * Deduces Zero<Real> from Zero()
   */
  Zero() -> Zero<Real>;

  /**
   * @brief Vector zero function
   * @tparam Scalar Scalar type for vector components
   * @ingroup ZeroSpecializations
   *
   * Represents the constant zero vector function:
   * @f[
   *   \mathbf{f}(x) = \mathbf{0} \quad \forall x \in \Omega
   * @f]
   * where @f$ \mathbf{0} \in \mathbb{R}^d @f$ is the zero vector in @f$ d @f$ dimensions.
   *
   * @see Zero
   */
  template <class Scalar>
  class Zero<Math::Vector<Scalar>> final
    : public VectorFunctionBase<Scalar, Zero<Math::Vector<Scalar>>>
  {
    public:
      /// @brief Vector value type
      using VectorType = Math::Vector<Scalar>;

      /// @brief Parent class type
      using Parent = VectorFunctionBase<Scalar, Zero<VectorType>>;

      /**
       * @brief Constructor with dimension
       * @param d Dimension of the zero vector
       */
      Zero(size_t d)
        : m_d(d)
      {}

      /**
       * @brief Copy constructor
       * @param other Zero vector function to copy from
       */
      Zero(const Zero& other)
        : Parent(other),
          m_d(other.m_d)
      {}

      /**
       * @brief Move constructor
       * @param other Zero vector function to move from
       */
      Zero(Zero&& other)
        : Parent(std::move(other)),
          m_d(std::move(other.m_d))
      {}

      /**
       * @brief Evaluates the zero vector function at a point
       * @param[in] p Point at which to evaluate (unused)
       * @returns Zero vector of dimension m_d
       */
      auto getValue(const Geometry::Point&) const
      {
        return VectorType::Zero(m_d);
      }

      /**
       * @brief Creates a copy of this zero vector function
       * @returns Pointer to newly allocated copy
       */
      Zero* copy() const noexcept override
      {
        return new Zero(*this);
      }

    private:
      const size_t m_d; ///< Dimension of the zero vector
  };

  /**
   * @brief CTAD for vector Zero
   *
   * Deduces Zero<Math::Vector<Real>> from Zero(size_t)
   */
  Zero(size_t) -> Zero<Math::Vector<Real>>;

  /**
   * @brief Convenience typedef for vector zero function
   * @tparam Scalar Scalar type for vector components
   */
  template <class Scalar>
  using VectorZero = Zero<Math::Vector<Scalar>>;
}

#endif

