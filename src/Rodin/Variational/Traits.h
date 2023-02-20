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
#include "Rodin/Utility/IsSpecialization.h"

#include "ForwardDecls.h"
#include "RangeType.h"

namespace Rodin::FormLanguage
{
  namespace Internal
  {
    template <bool B, bool S, bool V, bool M>
    struct RangeOfSAT
    {};

    template <>
    struct RangeOfSAT<true, false, false, false>
    {
      using Type = Boolean;
      Variational::RangeType Value = Variational::RangeType::Boolean;
    };

    template <>
    struct RangeOfSAT<false, true, false, false>
    {
      using Type = Scalar;
      Variational::RangeType Value = Variational::RangeType::Scalar;
    };

    template <>
    struct RangeOfSAT<false, false, true, false>
    {
      using Type = Math::Vector;
      Variational::RangeType Value = Variational::RangeType::Vector;
    };

    template <>
    struct RangeOfSAT<false, false, false, true>
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

  template <class T>
  struct IsBooleanRange
  {
    static constexpr const bool Value = std::is_same_v<T, Boolean>;
  };

  template <class T>
  struct IsScalarRange
  {
    static constexpr const bool Value = std::is_convertible_v<T, Scalar>;
  };

  template <class T>
  struct IsVectorRange
  {
    static constexpr const bool Value = std::is_assignable_v<T, Math::Vector>;
  };

  template <class T>
  struct IsMatrixRange
  {
    static constexpr const bool Value = std::is_assignable_v<T, Math::Matrix>;
  };

  template <class T>
  struct ResultOf
  {};

  template <class Derived>
  struct ResultOf<Variational::FunctionBase<Derived>>
  {
    using Type =
      std::invoke_result_t<Variational::FunctionBase<Derived>,
        const Geometry::Point&>;
  };

  template <class Derived, Variational::ShapeFunctionSpaceType Space>
  struct ResultOf<Variational::ShapeFunctionBase<Derived, Space>>
  {
    using Type =
      std::invoke_result_t<Variational::ShapeFunctionBase<Derived, Space>,
        Variational::ShapeComputator&, const Geometry::Point&>;
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
                                     IsScalarRange<ResultType>::Value,
                                     IsVectorRange<ResultType>::Value,
                                     IsMatrixRange<ResultType>::Value>
                                     ::Type;
  };

  template <class Derived, Variational::ShapeFunctionSpaceType Space>
  struct RangeOf<Variational::ShapeFunctionBase<Derived, Space>>
  {
    using ResultType = typename ResultOf<Variational::ShapeFunctionBase<Derived, Space>>::Type;
    static_assert(Utility::IsSpecialization<ResultType, Variational::TensorBasis>::Value);
    static_assert(HasValueTypeMember<ResultType>::Value);
    using ValueType = typename ResultType::ValueType;

    using Type =
      typename Internal::RangeOfSAT<IsBooleanRange<ValueType>::Value,
                                     IsScalarRange<ValueType>::Value,
                                     IsVectorRange<ValueType>::Value,
                                     IsMatrixRange<ValueType>::Value>
                                     ::Type;
  };

  template <class Derived, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::ShapeFunctionBase<Derived, Space>>
  {
    using ResultType = typename ResultOf<Variational::ShapeFunctionBase<Derived, Space>>::Type;
    using RangeType = typename RangeOf<Variational::ShapeFunctionBase<Derived, Space>>::Type;

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
