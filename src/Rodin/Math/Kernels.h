/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_KERNELS_H
#define RODIN_MATH_KERNELS_H

#include "Rodin/Math.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Math::Kernels
{
  template <class MatrixScalar, class VectorScalar, class DOFScalar>
  static void eliminate(
      SparseMatrix<MatrixScalar>& stiffness, Vector<VectorScalar>& mass,
      const IndexMap<DOFScalar>& dofs, size_t offset = 0)
  {
    auto* const valuePtr = stiffness.valuePtr();
    auto* const outerPtr = stiffness.outerIndexPtr();
    auto* const innerPtr = stiffness.innerIndexPtr();
    // Move essential degrees of freedom in the LHS to the RHS
    for (const auto& kv : dofs)
    {
      const Index& global = kv.first + offset;
      const auto& dof = kv.second;
      for (typename SparseMatrix<MatrixScalar>::InnerIterator it(stiffness, global); it; ++it)
         mass.coeffRef(it.row()) -= it.value() * dof;
    }
    for (const auto& [global, dof] : dofs)
    {
      // Impose essential degrees of freedom on RHS
      mass.coeffRef(global + offset) = dof;

      // Impose essential degrees of freedom on LHS
      for (auto i = outerPtr[global + offset]; i < outerPtr[global + offset + 1]; ++i)
      {
        assert(innerPtr[i] >= 0);
        // Assumes CCS format
        const Index row = innerPtr[i];
        valuePtr[i] = (row == global + offset);
        if (row != global + offset)
        {
          for (auto k = outerPtr[row]; k < outerPtr[row + 1]; k++)
          {
            if (static_cast<Index>(innerPtr[k]) == global + offset)
            {
               valuePtr[k] = 0;
               break;
            }
          }
        }
      }
    }
  }

  template <class Scalar>
  static void replace(
      const Vector<Scalar>& row,
      SparseMatrix<Scalar>& stiffness, Vector<Scalar>& mass,
      const IndexMap<Scalar>& dofs, size_t offset = 0)
  {
    for (const auto& kv : dofs)
    {
      const Index& global = kv.first + offset;
      const auto& dof = kv.second;
      mass.coeffRef(global) = dof;
      assert(row.size() >= 0);
      for (size_t i = 0; i < static_cast<size_t>(row.size()); i++)
        stiffness.insert(global, i) = row(i);
    }
  }
}

#endif

