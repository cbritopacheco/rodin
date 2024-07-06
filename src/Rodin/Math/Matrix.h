/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_MATRIX_H
#define RODIN_MATH_MATRIX_H

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Rodin/Configure.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Types.h"
#include "Rodin/Array.h"

namespace Rodin::Math
{
  /**
   * @brief Dense scalar valued matrix type.
   */
  template <class ScalarType>
  using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

  using RealMatrix = Matrix<Real>;

  /**
   * @brief Dense complex valued matrix type.
   */
  using ComplexMatrix = Matrix<Complex>;

  /**
   * @brief Spatial matrix
   */
  template <class ScalarType>
  using SpatialMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, 0,
        RODIN_MAXIMAL_SPACE_DIMENSION, RODIN_MAXIMAL_SPACE_DIMENSION>;

  /**
   * @brief Point matrix
   */
  using PointMatrix =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, 0, RODIN_MAXIMAL_SPACE_DIMENSION, Eigen::Dynamic>;

  /**
   * @brief Represents a fixed size matrix
   */
  template <class ScalarType, size_t Rows, size_t Cols>
  using FixedSizeMatrix = Eigen::Matrix<ScalarType, Rows, Cols>;
}

namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::Matrix<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number>
  struct Traits<Math::SpatialMatrix<Number>>
  {
    using ScalarType = Number;
  };
}

#endif

