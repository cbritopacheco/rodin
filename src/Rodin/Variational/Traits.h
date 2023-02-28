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
#include "Rodin/Utility/IsDetected.h"
#include "Rodin/Utility/IsSpecialization.h"

#include "ForwardDecls.h"
#include "RangeType.h"

namespace Rodin::FormLanguage
{
  namespace Internal
  {
    template <bool B, bool I, bool S, bool V, bool M>
    struct RangeOfSAT
    {};

    template <>
    struct RangeOfSAT<true, false, false, false, false>
    {
      using Type = Boolean;
      Variational::RangeType Value = Variational::RangeType::Boolean;
    };

    template <>
    struct RangeOfSAT<false, true, false, false, false>
    {
      using Type = Integer;
      Variational::RangeType Value = Variational::RangeType::Integer;
    };

    template <>
    struct RangeOfSAT<false, false, true, false, false>
    {
      using Type = Scalar;
      Variational::RangeType Value = Variational::RangeType::Scalar;
    };

    template <>
    struct RangeOfSAT<false, false, false, true, false>
    {
      using Type = Math::Vector;
      Variational::RangeType Value = Variational::RangeType::Vector;
    };

    template <>
    struct RangeOfSAT<false, false, false, false, true>
    {
      using Type = Math::Matrix;
      Variational::RangeType Value = Variational::RangeType::Matrix;
    };
  }

  template <typename T, typename = void>
  struct HasValueTypeMember
  {
    static constexpr const bool Value = false;
  };

  template <class T>
  struct HasValueTypeMember<T, std::void_t<typename T::ValueType>>
  {
    static constexpr const bool Value = true;
  };

  template <typename T, typename = void, typename = void>
  struct IsVectorAtCompileTime
  {
    static constexpr const bool Value = false;
  };

  template <class T>
  struct IsVectorAtCompileTime<T, std::void_t<typename T::IsVectorAtCompileTime>>
  {
    static constexpr const bool Value = T::IsVectorAtCompileTime;
  };

  template <class T>
  struct IsBooleanRange
  {
    static constexpr const bool Value = std::is_same_v<T, Boolean>;
  };

  template <class T>
  struct IsIntegerRange
  {
    static constexpr const bool Value = std::is_same_v<T, Integer>;
  };

  template <class T>
  struct IsScalarRange
  {
    static constexpr const bool Value = std::is_convertible_v<T, Scalar>;
  };

  template <class T, typename = void, typename = void>
  struct IsVectorRange
  {
    static constexpr const bool Value = false;
  };

  template <class T>
  struct IsVectorRange<T, std::void_t<decltype(T::RowsAtCompileTime)>, std::void_t<decltype(T::ColsAtCompileTime)>>
  {
    static constexpr const bool Value = T::ColsAtCompileTime == 1;
  };

  template <class T, typename = void, typename = void>
  struct IsMatrixRange
  {
    static constexpr const bool Value = false;
  };

  template <class T>
  struct IsMatrixRange<T, std::void_t<decltype(T::RowsAtCompileTime)>, std::void_t<decltype(T::ColsAtCompileTime)>>
  {
    static constexpr const bool Value =
      std::is_assignable_v<Math::Matrix, T>
      && !IsVectorRange<T,
        std::void_t<decltype(T::RowsAtCompileTime)>, std::void_t<decltype(T::ColsAtCompileTime)>>::Value;
  };

  template <class T>
  struct IsPlainObject;

  template <class EigenDerived>
  struct IsPlainObject<Eigen::MatrixBase<EigenDerived>>
  {
    static constexpr const bool Value =
      std::is_base_of_v<Eigen::PlainObjectBase<std::decay_t<EigenDerived>>, std::decay_t<EigenDerived>>;
  };

  template <class T>
  struct ResultOf
  {};

  template <class Derived>
  struct ResultOf<Variational::FunctionBase<Derived>>
  {
    using Type =
      std::invoke_result_t<Variational::FunctionBase<Derived>, const Geometry::Point&>;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using Type =
      std::invoke_result_t<Variational::ShapeFunctionBase<Derived, FES, Space>, const Geometry::Point&>;
  };

  template <class T>
  struct RangeOf
  {};

  template <class Derived>
  struct RangeOf<Variational::FunctionBase<Derived>>
  {
    using ResultType = typename ResultOf<Variational::FunctionBase<Derived>>::Type;

    using Type =
      typename Internal::RangeOfSAT<IsBooleanRange<ResultType>::Value,
                                    IsIntegerRange<ResultType>::Value,
                                    IsScalarRange<ResultType>::Value,
                                    IsVectorRange<ResultType>::Value,
                                    IsMatrixRange<ResultType>::Value>
                                    ::Type;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct RangeOf<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using ResultType = typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>::Type;
    static_assert(Utility::IsSpecialization<ResultType, Variational::TensorBasis>::Value,
        "EXPECTED TensorBasis AS RETURN TYPE. GOT SOMETHING ELSE. CHECK INSTANTIATION OF THE TEMPLATE CLASS.");
    static_assert(HasValueTypeMember<ResultType>::Value,
        "EXPECTED ValueType MEMBER DEFINITION DOES NOT EXIST.");
    using ValueType = typename ResultType::ValueType;

    using Type =
      typename Internal::RangeOfSAT<IsBooleanRange<ValueType>::Value,
                                    IsIntegerRange<ValueType>::Value,
                                    IsScalarRange<ValueType>::Value,
                                    IsVectorRange<ValueType>::Value,
                                    IsMatrixRange<ValueType>::Value>
                                    ::Type;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using ResultType = typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, Space>>::Type;
    using RangeType = typename RangeOf<Variational::ShapeFunctionBase<Derived, FES, Space>>::Type;

    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;
  };

  template <class Derived>
  struct Traits<Variational::FunctionBase<Derived>>
  {
    using ResultType = typename ResultOf<Variational::FunctionBase<Derived>>::Type;
    using RangeType = typename RangeOf<Variational::FunctionBase<Derived>>::Type;
  };
}

#endif
