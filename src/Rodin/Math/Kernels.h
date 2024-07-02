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
  inline
  static void eliminate(
      SparseMatrix<Real>& stiffness, Vector<Real>& mass,
      const IndexMap<Real>& dofs, size_t offset = 0)
  {
    Real* const valuePtr = stiffness.valuePtr();
    SparseMatrix<Real>::StorageIndex* const outerPtr = stiffness.outerIndexPtr();
    SparseMatrix<Real>::StorageIndex* const innerPtr = stiffness.innerIndexPtr();
    // Move essential degrees of freedom in the LHS to the RHS
    for (const auto& kv : dofs)
    {
      const Index& global = kv.first + offset;
      const auto& dof = kv.second;
      for (SparseMatrix<Real>::InnerIterator it(stiffness, global); it; ++it)
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
        valuePtr[i] = Real(row == global + offset);
        if (row != global + offset)
        {
          for (auto k = outerPtr[row]; 1; k++)
          {
            if (static_cast<Index>(innerPtr[k]) == global + offset)
            {
               valuePtr[k] = 0.0;
               break;
            }
          }
        }
      }
    }
  }

  inline
  static void replace(const Vector<Real>& row, SparseMatrix<Real>& stiffness, Vector<Real>& mass,
      const IndexMap<Real>& dofs, size_t offset = 0)
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

