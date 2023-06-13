/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_TENSOR_H
#define RODIN_MATH_TENSOR_H

#include <unsupported/Eigen/CXX11/Tensor>

#include "Rodin/Types.h"
#include "ForwardDecls.h"

namespace Rodin::Math
{
  /**
   * @brief Dense tensor type.
   * @tparam Rank Rank of tensor.
   */
  template <size_t Rank>
  using Tensor = Eigen::Tensor<Scalar, Rank>;

  /**
   * @brief Gets the tank of a tensor.
   */
  template <auto Rank>
  constexpr
  auto rank(const Tensor<Rank>& tensor)
  {
    return Rank;
  }

  template <class T, size_t Dim>
  class Slice;

  template <size_t Dim>
  class Slice<const Tensor<3>, Dim>
  {};

  template <>
  class Slice<const Tensor<3>, 0>
    : public Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  {
    public:
      using Parent =
        Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

      using Parent::operator=;

      Slice(const Tensor<3>& tensor, size_t offset)
        : Parent(tensor.data() + offset, tensor.dimension(1), tensor.dimension(2),
            { tensor.dimension(0) * tensor.dimension(1), tensor.dimension(0) })
      {}

      Slice(const Slice& other)
        : Parent(other)
      {}

      Slice(Slice&& other)
        : Parent(std::move(other))
      {}
  };

  template <>
  class Slice<const Tensor<3>, 1>
    : public Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  {
    public:
      using Parent =
        Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

      using Parent::operator=;

      Slice(const Tensor<3>& tensor, size_t offset)
        : Parent(tensor.data() + offset * tensor.dimension(0), tensor.dimension(0), tensor.dimension(2),
            { tensor.dimension(0) * tensor.dimension(1), 1 })
      {}

      Slice(const Slice& other)
        : Parent(other)
      {}

      Slice(Slice&& other)
        : Parent(std::move(other))
      {}
  };

  template <>
  class Slice<const Tensor<3>, 2>
    : public Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  {
    public:
      using Parent =
        Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

      using Parent::operator=;

      Slice(const Tensor<3>& tensor, size_t offset)
        : Parent(tensor.data() + offset * tensor.dimension(0) * tensor.dimension(1),
            tensor.dimension(0), tensor.dimension(1), { tensor.dimension(0), 1 })
      {}

      Slice(const Slice& other)
        : Parent(other)
      {}

      Slice(Slice&& other)
        : Parent(std::move(other))
      {}
  };

  template <>
  class Slice<Tensor<3>, 0>
    : public Eigen::Map<Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  {
    public:
      using Parent =
        Eigen::Map<Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

      using Parent::operator=;

      Slice(Tensor<3>& tensor, size_t offset)
        : Parent(tensor.data() + offset, tensor.dimension(1), tensor.dimension(2),
            { tensor.dimension(0) * tensor.dimension(1), tensor.dimension(0) })
      {}

      Slice(const Slice& other)
        : Parent(other)
      {}

      Slice(Slice&& other)
        : Parent(std::move(other))
      {}
  };

  template <>
  class Slice<Tensor<3>, 1>
    : public Eigen::Map<Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  {
    public:
      using Parent =
        Eigen::Map<Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

      using Parent::operator=;

      Slice(Tensor<3>& tensor, size_t offset)
        : Parent(tensor.data() + offset * tensor.dimension(0), tensor.dimension(0), tensor.dimension(2),
            { tensor.dimension(0) * tensor.dimension(1), 1 })
      {}

      Slice(const Slice& other)
        : Parent(other)
      {}

      Slice(Slice&& other)
        : Parent(std::move(other))
      {}
  };

  template <>
  class Slice<Tensor<3>, 2>
    : public Eigen::Map<Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  {
    public:
      using Parent =
        Eigen::Map<Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

      using Parent::operator=;

      Slice(Tensor<3>& tensor, size_t offset)
        : Parent(tensor.data() + offset * tensor.dimension(0) * tensor.dimension(1),
            tensor.dimension(0), tensor.dimension(1), { tensor.dimension(0), 1 })
      {}

      Slice(const Slice& other)
        : Parent(other)
      {}

      Slice(Slice&& other)
        : Parent(std::move(other))
      {}
  };
}

#endif
