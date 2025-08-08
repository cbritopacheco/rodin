/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_LINEARSYSTEM_H
#define RODIN_MATH_LINEARSYSTEM_H

#include "Rodin/Array.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class Operator, class Vector>
  struct Traits<Math::LinearSystem<Operator, Vector>>
  {
    using OperatorType = Operator;

    using VectorType = Vector;

    using ScalarType = typename Traits<OperatorType>::ScalarType;
  };
}

namespace Rodin::Math
{
  template <class Matrix, class Vector, class Derived>
  class LinearSystemBase
  {
    public:
      using MatrixType =
        Matrix;

      using VectorType =
        Vector;

      constexpr
      LinearSystemBase() = default;

      constexpr
      LinearSystemBase(const LinearSystemBase& other) = default;

      constexpr
      LinearSystemBase(LinearSystemBase&& other) noexcept = default;

      virtual ~LinearSystemBase() = default;

      constexpr
      LinearSystemBase& operator=(const LinearSystemBase& other) = default;

      constexpr
      LinearSystemBase& operator=(LinearSystemBase&& other) noexcept = default;

      template <class DOFScalar>
      constexpr
      LinearSystemBase& eliminate(
          const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        return static_cast<Derived&>(*this).eliminate(dofs, offset);
      }

      template <class DOFScalar>
      constexpr
      LinearSystemBase& merge(
          const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        return static_cast<Derived&>(*this).merge(dofs, offset);
      }

      constexpr
      MatrixType& getOperator()
      {
        return static_cast<Derived&>(*this).getOperator();
      }

      constexpr
      const MatrixType& getOperator() const
      {
        return static_cast<const Derived&>(*this).getOperator();
      }

      constexpr
      VectorType& getVector()
      {
        return static_cast<Derived&>(*this).getVector();
      }

      constexpr
      const VectorType& getVector() const
      {
        return static_cast<const Derived&>(*this).getVector();
      }

      constexpr
      VectorType& getSolution()
      {
        return static_cast<Derived&>(*this).getSolution();
      }

      constexpr
      const VectorType& getSolution() const
      {
        return static_cast<const Derived&>(*this).getSolution();
      }
  };

  template <class Matrix, class Vector>
  class LinearSystem;

  template <class MatrixScalar, class VectorScalar>
  class LinearSystem<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>>
    : public LinearSystemBase<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>, LinearSystem<MatrixScalar, VectorScalar>>
  {
    public:
      using MatrixType =
        Math::SparseMatrix<MatrixScalar>;

      using VectorType =
        Math::Vector<VectorScalar>;

      using Parent =
        LinearSystemBase<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>, LinearSystem<MatrixScalar, VectorScalar>>;

      constexpr
      LinearSystem() = default;

      constexpr
      LinearSystem(const LinearSystem& other)
        : Parent(other),
          m_operator(other.m_operator),
          m_vector(other.m_vector),
          m_solution(other.m_solution)
      {}

      constexpr
      LinearSystem(LinearSystem&& other) noexcept
        : Parent(std::move(other)),
          m_operator(std::move(other.m_operator)),
          m_vector(std::move(other.m_vector)),
          m_solution(std::move(other.m_solution))
      {}

      constexpr
      LinearSystem& operator=(const LinearSystem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_operator = other.m_operator;
          m_vector = other.m_vector;
          m_solution = other.m_solution;
        }
        return *this;
      }

      constexpr
      LinearSystem& operator=(LinearSystem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_operator = std::move(other.m_operator);
          m_vector = std::move(other.m_vector);
          m_solution = std::move(other.m_solution);
        }
        return *this;
      }

      virtual ~LinearSystem() = default;

      template <class DOFScalar>
      constexpr
      LinearSystem& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = this->getOperator();
        auto& mass = this->getVector();

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

      template <class DOFScalar>
      LinearSystem& merge(const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        auto& stiffness = this->getOperator();
        auto& mass = this->getVector();

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

      constexpr
      MatrixType& getOperator()
      {
        return m_operator;
      }

      constexpr
      const MatrixType& getOperator() const
      {
        return m_operator;
      }

      constexpr
      VectorType& getVector()
      {
        return m_vector;
      }

      constexpr
      const VectorType& getVector() const
      {
        return m_vector;
      }

      constexpr
      VectorType& getSolution()
      {
        return m_solution;
      }

      constexpr
      const VectorType& getSolution() const
      {
        return m_solution;
      }

    private:
      MatrixType m_operator; ///< The operator of the linear system.
      VectorType m_vector;   ///< The vector of the linear system.
      VectorType m_solution; ///< The solution vector of the linear system.
  };

  template <class MatrixScalar, class VectorScalar>
  class LinearSystem<Math::Matrix<MatrixScalar>, Math::Vector<VectorScalar>>
    : public LinearSystemBase<Math::Matrix<MatrixScalar>, Math::Vector<VectorScalar>, LinearSystem<MatrixScalar, VectorScalar>>
  {
    public:
      using MatrixType =
        Math::Matrix<MatrixScalar>;

      using VectorType =
        Math::Vector<VectorScalar>;

      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixScalar, VectorScalar>>;

      constexpr
      LinearSystem() = default;

      constexpr
      LinearSystem(const LinearSystem& other)
        : Parent(other),
          m_operator(other.m_operator),
          m_vector(other.m_vector),
          m_solution(other.m_solution)
      {}

      constexpr
      LinearSystem(LinearSystem&& other) noexcept
        : Parent(std::move(other)),
          m_operator(std::move(other.m_operator)),
          m_vector(std::move(other.m_vector)),
          m_solution(std::move(other.m_solution))
      {}

      virtual ~LinearSystem() = default;

      constexpr
      LinearSystem& operator=(const LinearSystem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_operator = other.m_operator;
          m_vector = other.m_vector;
          m_solution = other.m_solution;
        }
        return *this;
      }

      constexpr
      LinearSystem& operator=(LinearSystem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_operator = std::move(other.m_operator);
          m_vector = std::move(other.m_vector);
          m_solution = std::move(other.m_solution);
        }
        return *this;
      }

      template <class DOFScalar>
      LinearSystem& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = this->getOperator();
        auto& mass = this->getVector();

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

      template <class DOFScalar>
      LinearSystem& merge(const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        auto& stiffness = this->getOperator();
        auto& mass = this->getVector();

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

      constexpr
      MatrixType& getOperator()
      {
        return m_operator;
      }

      constexpr
      const MatrixType& getOperator() const
      {
        return m_operator;
      }

      constexpr
      VectorType& getVector()
      {
        return m_vector;
      }

      constexpr
      const VectorType& getVector() const
      {
        return m_vector;
      }

      constexpr
      VectorType& getSolution()
      {
        return m_solution;
      }

      constexpr
      const VectorType& getSolution() const
      {
        return m_solution;
      }

    private:
      MatrixType m_operator; ///< The operator of the linear system.
      VectorType m_vector;   ///< The vector of the linear system.
      VectorType m_solution; ///< The solution vector of the linear system.
  };
}

namespace std
{
  template <class Matrix, class Vector>
  struct tuple_size<Rodin::Math::LinearSystem<Matrix, Vector>>
    : std::integral_constant<std::size_t, 3> {};

  template <class Matrix, class Vector>
  struct tuple_element<0, Rodin::Math::LinearSystem<Matrix, Vector>>
  {
    using type = Matrix;
  };

  template <class Matrix, class Vector>
  struct tuple_element<1, Rodin::Math::LinearSystem<Matrix, Vector>>
  {
    using type = Vector;
  };

  template <class Matrix, class Vector>
  struct tuple_element<2, Rodin::Math::LinearSystem<Matrix, Vector>>
  {
    using type = Vector;
  };
}

namespace Rodin::Math
{
  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(LinearSystem<Matrix, Vector>& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0) return (ls.getOperator());
    else if constexpr (I == 1) return (ls.getSolution());
    else return (ls.getVector());
  }

  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(const LinearSystem<Matrix, Vector>& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0) return (ls.getOperator());
    else if constexpr (I == 1) return (ls.getSolution());
    else return (ls.getVector());
  }

  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(LinearSystem<Matrix, Vector>&& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0) return std::move(ls.getOperator());
    else if constexpr (I == 1) return std::move(ls.getSolution());
    else return std::move(ls.getVector());
  }

  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(const LinearSystem<Matrix, Vector>&& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0) return std::move(ls.getOperator());
    else if constexpr (I == 1) return std::move(ls.getSolution());
    else return std::move(ls.getVector());
  }
}

#endif

