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
        ::Mat& A = this->getOperator();
        ::Vec& b = this->getVector();

        std::vector<PetscInt> ids;
        ids.reserve(dofs.size());
        for (auto const& kv : dofs)
          ids.push_back(PetscInt(kv.first) + PetscInt(offset));
        IS is;
        ISCreateGeneral(PETSC_COMM_WORLD,
                        PetscInt(ids.size()), ids.data(),
                        PETSC_COPY_VALUES, &is);

        Vec x;
        VecDuplicate(b, &x);
        VecZeroEntries(x);
        for (auto const& kv : dofs)
        {
          const PetscInt  i   = PetscInt(kv.first) + PetscInt(offset);
          const PetscReal ui  = PetscReal(kv.second);
          VecSetValue(x, i, ui, INSERT_VALUES);
        }

        VecAssemblyBegin(x);
        VecAssemblyEnd(x);
        MatZeroRowsColumnsIS(A, is, 1.0, x, b);

        ISDestroy(&is);
        VecDestroy(&x);

        return *this;
      }

      template <class DOFScalar>
      LinearSystemBase& merge(
          const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        throw "Unimplemented.";

        return *this;
      }
  };
}

#endif
