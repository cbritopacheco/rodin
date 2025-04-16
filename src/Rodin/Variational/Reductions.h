/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLE_REDUCTIONS_H
#define RODIN_ASSEMBLE_REDUCTIONS_H

#include "Rodin/Math.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Variational::Reductions
{
  template <class MatrixScalar, class VectorScalar, class DOFScalar>
  static void eliminate(
      Math::SparseMatrix<MatrixScalar>& stiffness, Math::Vector<VectorScalar>& mass,
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
      for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(stiffness, global); it; ++it)
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
      const Math::Vector<Scalar>& row,
      Math::SparseMatrix<Scalar>& stiffness, Math::Vector<Scalar>& mass,
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

  template <class MatrixScalar, class VectorScalar, class DOFScalar>
  static void merge(
      Math::SparseMatrix<MatrixScalar>& stiffness, Math::Vector<VectorScalar>& mass,
      const IndexMap<DOFScalar>& dofs, size_t offset = 0)
  {
    std::deque<Index> q;
    IndexSet dependents;
    dependents.reserve(dofs.size());
    for (const auto& [k, v] : dofs)
      dependents.insert(v.first.begin(), v.first.end());

    for (auto it = dofs.begin(); it != dofs.end(); ++it)
    {
      const Index k = it->first;
      if (!dependents.contains(k))
        q.push_front(k);
    }

    // Perform breadth-first traversal
    while (q.size() > 0)
    {
      const Index parent = q.back();
      assert(stiffness.rows() >= 0);
      assert(stiffness.cols() >= 0);
      assert(parent < static_cast<size_t>(m_stiffness.rows()));
      assert(parent < static_cast<size_t>(m_stiffness.cols()));
      q.pop_back();

      auto find = dofs.find(parent);
      if (find == dofs.end())
        continue;

      const auto& [children, coeffs] = find->second;
      assert(children.size() > 0);

      for (const auto& child : children)
        q.push_front(child);

      assert(children.size() == coeffs.size());
      const size_t count = children.size();

      // Eliminate the parent column, adding it to the child columns
      for (size_t i = 0; i < count; i++)
      {
        const MatrixScalar coeff = coeffs.coeff(i);
        const Index child = children.coeff(i);
        stiffness.col(child) += coeff * stiffness.col(parent);
      }

      // Assumes CCS format
      for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(stiffness, parent); it; ++it)
        it.valueRef() = 0;

      // Eliminate the parent row, adding it to the child rows
      IndexMap<MatrixScalar> parentLookup;
      std::vector<IndexMap<MatrixScalar>> childrenLookup(children.size());
      for (size_t col = 0; col < static_cast<size_t>(stiffness.cols()); col++)
      {
        Boolean parentFound = false;
        size_t childrenFound = 0;
        for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(stiffness, col); it; ++it)
        {
          if (parentFound && childrenFound == count)
          {
            break;
          }
          else
          {
            const Index row = it.row();
            if (row == parent)
            {
              parentLookup[col] = it.value();
              it.valueRef() = 0;
              parentFound = true;
            }
            else
            {
              for (size_t i = 0; i < count; i++)
              {
                const Index child = children.coeff(i);
                if (row == child)
                {
                  childrenLookup[i][col] = it.value();
                  childrenFound += 1;
                }
              }
            }
          }
        }
      }

      for (const auto& [col, value] : parentLookup)
      {
        for (size_t i = 0; i < count; i++)
        {
          const MatrixScalar coeff = coeffs.coeff(i);
          childrenLookup[i][col] += coeff * value;
        }
      }

      for (size_t i = 0; i < count; i++)
      {
        const Index child = children.coeff(i);
        for (const auto& [col, value] : childrenLookup[i])
          stiffness.coeffRef(child, col) = value;
      }

      // Eliminate the parent entry, adding it to the child entries
      for (size_t i = 0; i < count; i++)
      {
        const MatrixScalar coeff = coeffs.coeff(i);
        const Index child = children.coeff(i);
        mass.coeffRef(child) += coeff * mass.coeff(parent);
      }
      mass.coeffRef(parent) = 0;
    }

    for (const auto& [parent, node] : dofs)
    {
      stiffness.coeffRef(parent, parent) = 1.0;
      const auto& [children, coeffs] = node;
      assert(children.size() >= 0);
      for (size_t i = 0; i < static_cast<size_t>(children.size()); i++)
      {
        const MatrixScalar coeff = coeffs.coeff(i);
        const Index child = children.coeff(i);
        stiffness.coeffRef(parent, child) = -coeff;
      }
    }
  }
}

#endif



