/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_SYSTEM_H
#define RODIN_PETSC_MATH_SYSTEM_H

#include <petsc.h>
#include "Rodin/Math/System.h"

namespace Rodin::Math
{
  template <>
  class LinearSystem<::Mat, ::Vec>
    : public LinearSystemBase<::Mat, ::Vec, LinearSystem<::Mat, ::Vec>>
  {
    public:
      using MatrixType = ::Mat;

      using VectorType = ::Vec;

      using Parent = LinearSystemBase<::Mat, ::Vec, LinearSystem<::Mat, ::Vec>>;

      using Parent::Parent;

      template <class DOFScalar>
      LinearSystem& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        return *this;
      }
  };
}

#endif
