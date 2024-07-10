/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRAITS_H
#define RODIN_VARIATIONAL_TRAITS_H

#include <type_traits>

#include "Rodin/Types.h"

#include "Rodin/Math/ForwardDecls.h"

#include "Rodin/Geometry/ForwardDecls.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Utility/HasValueMember.h"
#include "Rodin/Utility/IsSpecialization.h"

#include "ForwardDecls.h"
#include "RangeType.h"

namespace Rodin::FormLanguage
{
  template <class T>
  struct ResultOf;

  template <class Derived>
  struct ResultOf<Variational::FunctionBase<Derived>>
  {
    using Type =
      std::invoke_result_t<Variational::FunctionBase<Derived>, const Geometry::Point&>;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using Type = std::invoke_result_t<Variational::ShapeFunctionBase<Derived, FES, Space>, size_t>;
  };

  template <class T>
  struct RangeOf;

  template <>
  struct RangeOf<Boolean>
  {
    using Type = Boolean;
  };

  template <>
  struct RangeOf<Integer>
  {
    using Type = Integer;
  };

  template <>
  struct RangeOf<Real>
  {
    using Type = Real;
  };

  template <>
  struct RangeOf<Complex>
  {
    using Type = Complex;
  };

  template <class Scalar, int Rows, int Options, int MaxRows, int MaxCols>
  struct RangeOf<Eigen::Matrix<Scalar, Rows, 1, Options, MaxRows, MaxCols>>
  {
    using Type = Math::Vector<Scalar>;
  };

  template <class Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  struct RangeOf<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>
  {
    using Type = Math::Matrix<Scalar>;
  };

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

  template <class Derived>
  struct RangeOf<Variational::FunctionBase<Derived>>
  {
    using ResultType = typename ResultOf<Variational::FunctionBase<Derived>>::Type;
    using Type = typename RangeOf<ResultType>::Type;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct RangeOf<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using ResultType = typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>::Type;
    using Type = typename RangeOf<ResultType>::Type;
  };
}

#endif
