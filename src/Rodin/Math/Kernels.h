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
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

namespace Rodin::Math::Kernels
{
  inline
  static void add(
      Math::Matrix& out, const Math::Matrix& in,
      const IndexArray& rows, const IndexArray& cols)
  {
    assert(rows.size() >= 0);
    assert(cols.size() >= 0);
    assert(in.rows() == rows.size());
    assert(in.cols() == cols.size());
    out(rows, cols).noalias() += in;
  }

  inline
  static void add(Math::Vector& out, const Math::Vector& in, const IndexArray& s)
  {
    assert(in.size() == s.size());
    out(s).noalias() += in;
  }

  inline
  static void add(
      std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
      const IndexArray& rows, const IndexArray& cols)
  {
    assert(rows.size() >= 0);
    assert(cols.size() >= 0);
    assert(in.rows() == rows.size());
    assert(in.cols() == cols.size());
    for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
    {
      for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
      {
        const Scalar s = in(i, j);
        if (s != Scalar(0))
          out.emplace_back(rows(i), cols(j), s);
      }
    }
  }

  inline
  static void eliminate(SparseMatrix& stiffness, Vector& mass,
      const IndexMap<Scalar>& dofs, size_t offset = 0)
  {
    Scalar* const valuePtr = stiffness.valuePtr();
    Math::SparseMatrix::StorageIndex* const outerPtr = stiffness.outerIndexPtr();
    Math::SparseMatrix::StorageIndex* const innerPtr = stiffness.innerIndexPtr();
    // Move essential degrees of freedom in the LHS to the RHS
    for (const auto& kv : dofs)
    {
      const Index& global = kv.first + offset;
      const auto& dof = kv.second;
      for (Math::SparseMatrix::InnerIterator it(stiffness, global); it; ++it)
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
        valuePtr[i] = Scalar(row == global + offset);
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
  static void replace(const Math::Vector& row, SparseMatrix& stiffness, Vector& mass,
      const IndexMap<Scalar>& dofs, size_t offset = 0)
  {
    for (const auto& kv : dofs)
    {
      const Index& global = kv.first + offset;
      const auto& dof = kv.second;
      mass.coeffRef(global) = dof;
      for (size_t i = 0; i < row.size(); i++)
        stiffness.insert(global, i) = row(i);
    }
  }
}

#endif

