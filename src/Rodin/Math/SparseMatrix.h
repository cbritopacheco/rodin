/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_SPARSEMATRIX_H
#define RODIN_MATH_SPARSEMATRIX_H

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Rodin::Math
{
  class SparseMatrix : public Eigen::SparseMatrix<double>
  {
   public:
    using Parent = Eigen::SparseMatrix<double>;

    SparseMatrix() = default;

    template <typename OtherDerived>
    SparseMatrix(const Eigen::MatrixBase<OtherDerived>& other)
      : Parent(other)
    {}

    template <typename OtherDerived>
    SparseMatrix& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
       this->Parent::operator=(other);
       return *this;
    }
  };
}

#endif

