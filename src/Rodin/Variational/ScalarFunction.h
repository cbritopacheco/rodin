/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ScalarFunction.h
 * @brief Base class for scalar-valued functions in variational formulations.
 *
 * This file defines ScalarFunctionBase, which extends FunctionBase to represent
 * functions mapping points in the domain to scalar values @f$ f: \Omega \to \mathbb{F} @f$,
 * where @f$ \mathbb{F} @f$ is a scalar field (typically Real or Complex).
 */
#ifndef RODIN_VARIATIONAL_SCALARFUNCTION_H
#define RODIN_VARIATIONAL_SCALARFUNCTION_H

#include "ForwardDecls.h"

#include "Function.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Derived>
  struct Traits<Variational::ScalarFunctionBase<Scalar, Derived>>
  {
    using ScalarType = Scalar;
    using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup RealFunctionSpecializations RealFunction Template Specializations
   * @brief Template specializations of the RealFunction class.
   * @see RealFunction
   */

  /**
   * @brief Base class for scalar-valued functions with templated scalar type.
   *
   * ScalarFunctionBase extends FunctionBase to represent functions mapping
   * points in the domain to scalar values:
   * @f[
   *    f: \Omega \to \mathbb{F}
   * @f]
   * where @f$ \mathbb{F} @f$ is a scalar field (typically Real, Complex, or Boolean).
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex, Boolean)
   * @tparam Derived The derived class following CRTP pattern
   *
   * @note This serves as an intermediate base between FunctionBase and concrete
   * scalar function types like RealFunction, ComplexFunction, and BooleanFunction.
   *
   * @see FunctionBase, RealFunction, ComplexFunction, BooleanFunction
   */
  template <class Scalar, class Derived>
  class ScalarFunctionBase
    : public FunctionBase<ScalarFunctionBase<Scalar, Derived>>
  {
    public:
      /// @brief Type of scalar values returned by this function
      using ScalarType = Scalar;

      /// @brief Parent class type
      using Parent = FunctionBase<ScalarFunctionBase<ScalarType, Derived>>;

      /// @brief Import traceOf methods from parent
      using Parent::traceOf;

      /// @brief Import operator() from parent
      using Parent::operator();

      /// @brief Default constructor
      ScalarFunctionBase() = default;

      /// @brief Copy constructor
      /// @param[in] other Function to copy from
      ScalarFunctionBase(const ScalarFunctionBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor
      /// @param[in] other Function to move from
      ScalarFunctionBase(ScalarFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /// @brief Virtual destructor
      virtual ~ScalarFunctionBase() = default;

      /**
       * @brief Evaluates the scalar function at a point.
       *
       * CRTP method that delegates to the derived class's implementation.
       *
       * @param[in] p Point at which to evaluate the function
       * @returns Scalar value at the given point
       */
      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Creates a polymorphic copy of the function.
       * @returns Pointer to a newly allocated copy
       * @note Caller is responsible for managing the returned pointer
       */
      virtual ScalarFunctionBase* copy() const noexcept override = 0;
  };
}

#endif

