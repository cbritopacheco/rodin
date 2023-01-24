/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_DENSEMATRIX_H
#define RODIN_MATH_DENSEMATRIX_H

#include <Eigen/Core>
#include <Eigen/Dense>

namespace Rodin::Math
{
  class DenseMatrix : public Eigen::MatrixX<double>
  {
   public:
    using Parent = Eigen::MatrixX<double>;

    DenseMatrix() = default;

    template <typename OtherDerived>
    DenseMatrix(const Eigen::MatrixBase<OtherDerived>& other)
      : Parent(other)
    {}

    template <typename OtherDerived>
    DenseMatrix& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
       this->Parent::operator=(other);
       return *this;
    }
  };
}

#endif
