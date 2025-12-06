/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRAITS_H
#define RODIN_VARIATIONAL_TRAITS_H

/**
 * @file
 * @brief Type traits for variational formulations
 *
 * This file defines type traits for extracting result types and range types
 * from variational objects like functions and shape functions. These traits
 * enable compile-time type deduction for form language expressions.
 */

#include <type_traits>

#include "Rodin/Types.h"

#include "Rodin/Math/ForwardDecls.h"

#include "Rodin/Geometry/ForwardDecls.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Utility/HasValueMember.h"
#include "Rodin/Utility/IsSpecialization.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type trait to deduce the result type of an expression
   * @tparam T Type to deduce result from
   *
   * Provides a member typedef `Type` that represents the return type
   * of evaluating the expression represented by `T`.
   */
  template <class T>
  struct ResultOf;

  /**
   * @brief Specialization of ResultOf for FunctionBase
   * @tparam Derived Derived function class
   *
   * Deduces the result type when evaluating a function at a geometric point.
   */
  template <class Derived>
  struct ResultOf<Variational::FunctionBase<Derived>>
  {
    using Type =
      std::invoke_result_t<Variational::FunctionBase<Derived>, const Geometry::Point&>;
  };

  /**
   * @brief Specialization of ResultOf for ShapeFunctionBase
   * @tparam Derived Derived shape function class
   * @tparam FES Finite element space type
   * @tparam Space Shape function space type (Trial or Test)
   *
   * Deduces the result type when evaluating a shape function with a degree of freedom index.
   */
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using Type =
      std::invoke_result_t<Variational::ShapeFunctionBase<Derived, FES, Space>, size_t>;
  };

  /**
   * @brief Type trait to deduce the mathematical range type of an expression
   * @tparam T Type to deduce range from
   *
   * Provides a member typedef `Type` that represents the mathematical range
   * (e.g., Real, Vector, Matrix) of the expression represented by `T`.
   */
  template <class T>
  struct RangeOf;

  /**
   * @brief Range type for Boolean values
   */
  template <>
  struct RangeOf<Boolean>
  {
    using Type =
      Boolean;
  };

  /**
   * @brief Range type for Integer values
   */
  template <>
  struct RangeOf<Integer>
  {
    using Type =
      Integer;
  };

  /**
   * @brief Range type for Real values
   */
  template <>
  struct RangeOf<Real>
  {
    using Type =
      Real;
  };

  /**
   * @brief Range type for Complex values
   */
  template <>
  struct RangeOf<Complex>
  {
    using Type =
      Complex;
  };

  /**
   * @brief Range type for Eigen column vectors
   * @tparam Scalar Scalar type
   * @tparam Rows Number of rows
   * @tparam Options Eigen storage options
   * @tparam MaxRows Maximum number of rows
   * @tparam MaxCols Maximum number of columns
   */
  template <class Scalar, int Rows, int Options, int MaxRows, int MaxCols>
  struct RangeOf<Eigen::Matrix<Scalar, Rows, 1, Options, MaxRows, MaxCols>>
  {
    using Type =
      Math::Vector<Scalar>;
  };

  /**
   * @brief Range type for Eigen matrices
   * @tparam Scalar Scalar type
   * @tparam Rows Number of rows
   * @tparam Cols Number of columns
   * @tparam Options Eigen storage options
   * @tparam MaxRows Maximum number of rows
   * @tparam MaxCols Maximum number of columns
   */
  template <class Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  struct RangeOf<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>
  {
    using Type =
      Math::Matrix<Scalar>;
  };

  /**
   * @brief General range type deduction for Eigen expression templates
   * @tparam MatrixXpr Eigen expression type
   *
   * Deduces whether the expression is a vector or matrix based on compile-time properties.
   */
  template <class MatrixXpr>
  struct RangeOf
  {
    using Type =
      std::conditional_t<
        MatrixXpr::IsVectorAtCompileTime, Math::Vector<typename MatrixXpr::Scalar>,
        std::conditional_t<
          MatrixXpr::ColsAtCompileTime == 1, Math::Vector<typename MatrixXpr::Scalar>,
          Math::Matrix<typename MatrixXpr::Scalar>
        >
      >;
  };

  /**
   * @brief Range type for FunctionBase
   * @tparam Derived Derived function class
   *
   * Deduces the range by first determining the result type, then mapping it to a range type.
   */
  template <class Derived>
  struct RangeOf<Variational::FunctionBase<Derived>>
  {
    using ResultType =
      typename ResultOf<Variational::FunctionBase<Derived>>::Type;
    using Type =
      typename RangeOf<std::remove_cvref_t<ResultType>>::Type;
  };

  /**
   * @brief Range type for ShapeFunctionBase
   * @tparam Derived Derived shape function class
   * @tparam FES Finite element space type
   * @tparam Space Shape function space type
   *
   * Deduces the range by first determining the result type, then mapping it to a range type.
   */
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct RangeOf<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using ResultType =
      typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>::Type;
    using Type =
      typename RangeOf<std::remove_cvref_t<ResultType>>::Type;
  };
}

#endif
