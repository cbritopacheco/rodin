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
  template <size_t Rank>
  using Tensor = Eigen::Tensor<Scalar, Rank>;

  template <auto Rank>
  constexpr
  auto rank(const Tensor<Rank>& tensor)
  {
    return Rank;
  }

  inline
  Eigen::Map<const Math::Matrix, Tensor<3>::Options, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  slice(const Tensor<3>& tensor, size_t dim, size_t offset)
  {
    if (dim == 0)
    {
      return {
        tensor.data() + offset,
        tensor.dimension(1), tensor.dimension(2),
        { tensor.dimension(0) * tensor.dimension(1), tensor.dimension(0) }
      };
    }
    else if (dim == 1)
    {
      return {
        tensor.data() + offset * tensor.dimension(0),
        tensor.dimension(0), tensor.dimension(2),
        { tensor.dimension(0) * tensor.dimension(1), 1 }
      };
    }
    else if (dim == 2)
    {
      return {
        tensor.data() + offset * tensor.dimension(0) * tensor.dimension(1),
        tensor.dimension(0), tensor.dimension(1),
        { tensor.dimension(0), 1 }
      };
    }
    else
    {
      assert(dim < 3);
    }
    assert(false);
    return { nullptr, 0, 0 };
  }

  inline
  Eigen::Map<Math::Matrix, Tensor<3>::Options, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  slice(Tensor<3>& tensor, size_t dim, size_t offset)
  {
    if (dim == 0)
    {
      return {
        tensor.data() + offset,
        tensor.dimension(1), tensor.dimension(2),
        { tensor.dimension(0) * tensor.dimension(1), tensor.dimension(0) }
      };
    }
    else if (dim == 1)
    {
      return {
        tensor.data() + offset * tensor.dimension(0),
        tensor.dimension(0), tensor.dimension(2),
        { tensor.dimension(0) * tensor.dimension(1), 1 }
      };
    }
    else if (dim == 2)
    {
      return {
        tensor.data() + offset * tensor.dimension(0) * tensor.dimension(1),
        tensor.dimension(0), tensor.dimension(1),
        { tensor.dimension(0), 1 }
      };
    }
    else
    {
      assert(dim < 3);
    }
    assert(false);
    return { nullptr, 0, 0 };
  }
}

#endif
