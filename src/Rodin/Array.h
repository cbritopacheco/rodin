/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ARRAY_H
#define RODIN_ARRAY_H

#include <Eigen/Core>
#include <Eigen/Dense>

namespace Rodin
{
  template <class Scalar>
  using Array = Eigen::ArrayX<Scalar>;
}

#endif

