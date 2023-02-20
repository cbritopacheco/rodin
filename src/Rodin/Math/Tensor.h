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

namespace Rodin::Math
{
  template <size_t Rank>
  using Tensor = Eigen::Tensor<Scalar, Rank>;
}

#endif
