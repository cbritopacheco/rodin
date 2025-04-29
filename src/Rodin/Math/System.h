/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_SYSTEM_H
#define RODIN_MATH_SYSTEM_H

#include "Rodin/Math.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Math
{
  template <class Matrix, class Vector>
  class System;

  template <class Matrix, class Vector>
  System(Matrix&, Vector&) -> System<Matrix, Vector>;

  template <class MatrixScalar, class VectorScalar>
  class System<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>>
  {
    public:
      using MatrixType = Math::SparseMatrix<MatrixScalar>;
      using VectorType = Math::Vector<VectorScalar>;

      System(MatrixType& stiffness, VectorType& mass)
        : m_stiffness(stiffness), m_mass(mass)
      {}

      template <class DOFScalar>
      System& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = m_stiffness.get();
        auto& mass = m_mass.get();

        auto* const valuePtr = stiffness.valuePtr();
        auto* const outerPtr = stiffness.outerIndexPtr();
        auto* const innerPtr = stiffness.innerIndexPtr();
        // Move essential degrees of freedom in the LHS to the RHS
        for (const auto& kv : dofs)
        {
          const Index& global = kv.first;
          const auto& dof = kv.second;
          for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(stiffness, global + offset); it; ++it)
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
        return *this;
      }

      template <class Row, class DOFScalar>
      System& replace(const Row& row, const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = m_stiffness.get();
        auto& mass = m_mass.get();
        for (const auto& kv : dofs)
        {
          const Index& global = kv.first;
          const auto& dof = kv.second;
          mass.coeffRef(global + offset) = dof;
          assert(row.size() >= 0);
          for (size_t i = 0; i < static_cast<size_t>(row.size()); i++)
            stiffness.insert(global + offset, i) = row(i);
        }
        return *this;
      }

      template <class DOFScalar>
      System& merge(const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        auto& stiffness = m_stiffness.get();
        auto& mass = m_mass.get();

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
          assert(parent < static_cast<size_t>(stiffness.rows()));
          assert(parent < static_cast<size_t>(stiffness.cols()));
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
        return *this;
      }

    private:
      std::reference_wrapper<MatrixType> m_stiffness;
      std::reference_wrapper<VectorType> m_mass;
  };

  template <class MatrixScalar, class VectorScalar>
  class System<Math::Matrix<MatrixScalar>, Math::Vector<VectorScalar>>
  {
    public:
      using MatrixType = Math::Matrix<MatrixScalar>;

      using VectorType = Math::Vector<VectorScalar>;

      System(MatrixType& stiffness, VectorType& mass)
        : m_stiffness(stiffness), m_mass(mass)
      {}

      template <class DOFScalar>
      System& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = m_stiffness.get();
        auto& mass = m_mass.get();

        // Move essential degrees of freedom in the LHS to the RHS
        for (const auto& kv : dofs)
        {
          const Index& global = kv.first;
          const auto& dof = kv.second;
          mass -= dof * stiffness.col(global + offset);
        }
        for (const auto& [global, dof] : dofs)
        {
          // Impose essential degrees of freedom on RHS
          mass.coeffRef(global + offset) = dof;

          // Impose essential degrees of freedom on LHS
          stiffness.col(global + offset).setZero();
          stiffness.row(global + offset).setZero();
          stiffness.coeffRef(global + offset, global + offset) = 1;
        }

        return *this;
      }

      /**
       * @brief  Replaces specified rows in a dense stiffness matrix and entries in a mass vector.
       *
       * This function writes Dirichlet‐type boundary conditions (or any prescribed values)
       * into a dense stiffness matrix and a mass vector.  For each degree of freedom (DoF)
       * listed in @p dofs, it sets the corresponding entry in @p mass to the prescribed
       * value, and overwrites the entire row in @p stiffness with the values provided in
       * @p row.
       *
       * @tparam Scalar     Numeric scalar type (e.g., double, float).
       *
       * @param[in,out] stiffness  Dense stiffness matrix to be modified.  Row @c (global+offset)
       *                            is overwritten with the contents of @p row.
       * @param[in,out] mass       Mass (or load) vector to be modified.  Entry
       *                            @c (global+offset) is set to the prescribed DoF value.
       * @param[in]     row        Vector of length N containing the new row values for
       *                            each column 0..N−1 in @p stiffness.
       * @param[in]     dofs       Mapping from global index to prescribed value:
       *                            - Key:   global row index (before @p offset)
       *                            - Value: the value to write into @p mass at that index
       * @param[in]     offset     Optional row‐index offset to apply to all global indices.
       *                            Defaults to 0.
       *
       * @pre  For every entry in @p dofs, @f$ 0 \leq \text{global} + \text{offset} < stiffness.rows() @f$.
       * @pre  row.size() == stiffness.cols().
       */
      template <class Row, class DOFScalar>
      System& replace(const Row& row, const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = m_stiffness.get();
        auto& mass = m_mass.get();
        for (const auto& kv : dofs)
        {
          const Index& global = kv.first;
          const auto& dof = kv.second;
          mass(global + offset) = dof;
          assert(row.size() > 0);
          for (size_t i = 0; i < static_cast<size_t>(row.size()); ++i)
            stiffness(global + offset, i) = row(i);
        }
        return *this;
      }

      template <class DOFScalar>
      System& merge(const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        auto& stiffness = m_stiffness.get();
        auto& mass = m_mass.get();

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
          assert(parent < static_cast<size_t>(stiffness.rows()));
          assert(parent < static_cast<size_t>(stiffness.cols()));
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
            const DOFScalar coeff = coeffs.coeff(i);
            const Index child = children.coeff(i);
            stiffness.col(child) += coeff * stiffness.col(parent);
          }

          stiffness.col(parent).setZero();

          // Eliminate the parent row, adding it to the child rows
          IndexMap<MatrixScalar> parentLookup;
          std::vector<IndexMap<MatrixScalar>> childrenLookup(children.size());
          for (size_t col = 0; col < static_cast<size_t>(stiffness.cols()); col++)
          {
            bool parentFound = false;
            size_t childrenFound = 0;
            for (typename Math::Matrix<MatrixScalar>::InnerIterator it(stiffness, col); it; ++it)
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
                  stiffness.coeffRef(it.row(), it.col()) = 0;
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
              const DOFScalar coeff = coeffs.coeff(i);
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
            const DOFScalar coeff = coeffs.coeff(i);
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
            const DOFScalar coeff = coeffs.coeff(i);
            const Index child = children.coeff(i);
            stiffness.coeffRef(parent, child) = -coeff;
          }
        }
        return *this;
      }

    private:
      std::reference_wrapper<MatrixType> m_stiffness;
      std::reference_wrapper<VectorType> m_mass;
  };
}

#endif



