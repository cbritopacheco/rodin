/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_LINEARSYSTEM_H
#define RODIN_MATH_LINEARSYSTEM_H

#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Array.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type traits for a Math::LinearSystem: exposes the operator,
   * vector, and scalar types.
   */
  template <class Operator, class Vector>
  struct Traits<Math::LinearSystem<Operator, Vector>>
  {
    /// @brief Assembled operator type.
    using OperatorType = Operator;

    /// @brief Vector type of the linear system.
    using VectorType = Vector;

    /// @brief Scalar value type.
    using ScalarType = typename Traits<OperatorType>::ScalarType;
  };
}

namespace Rodin::Math
{
  /**
   * @brief Base class for linear systems of the form @f$ Ax = b @f$.
   *
   * This class provides a CRTP (Curiously Recurring Template Pattern) base 
   * for linear systems where @f$ A @f$ is a matrix operator, @f$ x @f$ is 
   * the solution vector, and @f$ b @f$ is the right-hand side vector.
   *
   * @tparam Matrix Matrix type (e.g., SparseMatrix, Matrix)
   * @tparam Vector Vector type
   * @tparam Derived Derived class implementing the CRTP pattern
   */
  template <class Matrix, class Vector, class Derived>
  class LinearSystemBase
  {
    public:
      /// Matrix type for the linear operator
      using MatrixType = Matrix;

      /// Vector type for solution and right-hand side
      using VectorType = Vector;

      /// Default constructor
      constexpr LinearSystemBase() = default;

      /**
       * Copy constructor
       * @param other Object to copy from.
       */
      constexpr LinearSystemBase(const LinearSystemBase& other) = default;

      /**
       * Move constructor
       * @param other Object to move from.
       */
      constexpr LinearSystemBase(LinearSystemBase&& other) noexcept = default;

      /// Virtual destructor
      virtual ~LinearSystemBase() = default;

      /**
       * Copy assignment operator
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      constexpr LinearSystemBase& operator=(const LinearSystemBase& other) = default;

      /**
       * Move assignment operator
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      constexpr LinearSystemBase& operator=(LinearSystemBase&& other) noexcept = default;

      /**
       * @brief Eliminates degrees of freedom from the linear system.
       * @tparam DOFScalar Scalar type for DOF values
       * @param dofs Map of DOF indices to values to eliminate
       * @param offset Optional offset for DOF indices
       * @return Reference to this linear system
       */
      template <class DOFScalar>
      constexpr LinearSystemBase& eliminate(
        const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        return static_cast<Derived&>(*this).eliminate(dofs, offset);
      }

      /**
       * @brief Merges additional contributions into the linear system.
       * @tparam DOFScalar Scalar type for DOF values
       * @param dofs Map of DOF indices to index arrays and vector contributions
       * @param offset Optional offset for DOF indices
       * @return Reference to this linear system
       */
      template <class DOFScalar>
      constexpr LinearSystemBase& merge(
        const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs,
        size_t offset = 0)
      {
        return static_cast<Derived&>(*this).merge(dofs, offset);
      }

      /**
       * @brief Gets the matrix operator @f$ A @f$ of the linear system.
       * @return Reference to the matrix operator
       */
      /// @brief Returns the system operator matrix @f$ A @f$.
      constexpr MatrixType& getOperator()
      {
        return static_cast<Derived&>(*this).getOperator();
      }

      /**
       * @brief Gets the matrix operator @f$ A @f$ of the linear system.
       * @return Const reference to the matrix operator
       */
      /// @brief Returns the system operator matrix @f$ A @f$.
      constexpr const MatrixType& getOperator() const
      {
        return static_cast<const Derived&>(*this).getOperator();
      }

      /**
       * @brief Gets the right-hand side vector @f$ b @f$ of the linear system.
       * @return Reference to the right-hand side vector
       */
      /// @brief Returns the right-hand side vector @f$ b @f$.
      constexpr VectorType& getVector()
      {
        return static_cast<Derived&>(*this).getVector();
      }

      /**
       * @brief Gets the right-hand side vector @f$ b @f$ of the linear system.
       * @return Const reference to the right-hand side vector
       */
      /// @brief Returns the right-hand side vector @f$ b @f$.
      constexpr const VectorType& getVector() const
      {
        return static_cast<const Derived&>(*this).getVector();
      }

      /**
       * @brief Gets the solution vector @f$ x @f$ of the linear system.
       * @return Reference to the solution vector
       */
      /// @brief Returns the solution vector @f$ x @f$.
      constexpr VectorType& getSolution()
      {
        return static_cast<Derived&>(*this).getSolution();
      }

      /**
       * @brief Gets the solution vector @f$ x @f$ of the linear system.
       * @return Const reference to the solution vector
       */
      /// @brief Returns the solution vector @f$ x @f$.
      constexpr const VectorType& getSolution() const
      {
        return static_cast<const Derived&>(*this).getSolution();
      }
  };

  /**
   * @brief Generic linear system template declaration.
   * @tparam Matrix Matrix type
   * @tparam Vector Vector type
   */
  template <class Matrix, class Vector>
  class LinearSystem;

  /**
   * @brief Linear system specialization for sparse matrices.
   *
   * This specialization handles linear systems @f$ Ax = b @f$ where @f$ A @f$ is a 
   * sparse matrix. This is commonly used for large systems arising from finite 
   * element discretizations.
   *
   * @tparam MatrixScalar Scalar type for matrix elements
   * @tparam VectorScalar Scalar type for vector elements
   */
  template <class MatrixScalar, class VectorScalar>
  class LinearSystem<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>>
    : public LinearSystemBase<Math::SparseMatrix<MatrixScalar>,
        Math::Vector<VectorScalar>,
        LinearSystem<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>>>
  {
    public:
      /// Sparse matrix type for the linear operator
      using MatrixType = Math::SparseMatrix<MatrixScalar>;

      /// Vector type for solution and right-hand side
      using VectorType = Math::Vector<VectorScalar>;

      /// Parent class type
      using Parent =
        LinearSystemBase<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>,
          LinearSystem<Math::SparseMatrix<MatrixScalar>, Math::Vector<VectorScalar>>>;

      /// Default constructor
      constexpr LinearSystem() = default;

      /**
       * Copy constructor
       * @param other Object to copy from.
       */
      constexpr LinearSystem(const LinearSystem& other)
        : Parent(other),
          m_operator(other.m_operator),
          m_vector(other.m_vector),
          m_solution(other.m_solution)
      {}

      /**
       * Move constructor
       * @param other Object to move from.
       */
      constexpr LinearSystem(LinearSystem&& other) noexcept
        : Parent(std::move(other)),
          m_operator(std::move(other.m_operator)),
          m_vector(std::move(other.m_vector)),
          m_solution(std::move(other.m_solution))
      {}

      /**
       * Copy assignment operator
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      constexpr LinearSystem& operator=(const LinearSystem& other)
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

      /**
       * Move assignment operator
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      constexpr LinearSystem& operator=(LinearSystem&& other) noexcept
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

      /// Virtual destructor
      virtual ~LinearSystem() = default;

      /**
       * @brief Imposes prescribed coefficients on a square sparse system.
       *
       * # Architecture
       *
       * Let @f$C@f$ be the constrained indices, @f$Q@f$ their diagonal
       * indicator, @f$P=I-Q@f$, and @f$g@f$ the prescribed vector extended
       * by zero outside @f$C@f$. The operation replaces @f$Ax=b@f$ by
       * @f[
       *   \widetilde A=PAP+Q,\qquad
       *   \widetilde b=P(b-Ag)+g.
       * @f]
       * Thus free equations retain their original coefficients and receive
       * the lifting @f$-A_{FC}g_C@f$; constrained equations become
       * @f$x_C=g_C@f$. Neither symmetry of values nor transpose storage is
       * required. In particular, no conjugation enters the lifting for complex
       * systems. The solution vector is preserved for warm starts.
       *
       * # Implementation
       *
       * All constrained columns contribute to the load before any coefficient
       * is cleared. Their diagonal positions are recorded after each column's
       * traversal. Absent diagonals are inserted together by one backwards
       * merge of the compressed arrays, preserving sorted row indices and all
       * existing stored entries, including explicit zeros.
       *
       * A byte per row identifies constraints. Adjacent constrained columns
       * are cleared as contiguous storage ranges. Between these ranges,
       * blocks of sixteen row indices form a bit mask: empty masks require no
       * coefficient stores, while nonempty masks visit only their set bits.
       * Sixteen is a kernel granularity, not a numerical threshold. The tail
       * uses bit masking for Real and conditional zero assignment for other
       * scalars; multiplication by zero would not clear nonfinite values.
       * Finally the diagonals and prescribed load entries are restored.
       *
       * For a single contiguous constraint interval, a bounded sample may
       * select exact row-range pruning. This heuristic requires at least eight
       * entries per column on average and at least half of sixteen storage-
       * quantile samples to reject the interval. These measured policy choices
       * select a kernel only: each skipped column is independently certified
       * by its actual row bounds. Partial columns use binary searches.
       * Fragmented constraints always use the block kernel.
       *
       * With @f$n@f$ rows, @f$m@f$ stored entries and @f$k@f$ constraints,
       * the block path takes @f$O(n+m+k)@f$ work and @f$O(n+k)@f$ scratch.
       * The optional pruning path additionally performs at most two binary
       * searches per overlapping free column. Matrix storage grows only by
       * the number of absent constrained diagonals. No transpose index, full
       * matrix copy, persistent cache, or worker pool is constructed.
       *
       * @tparam DOFScalar Prescribed scalar type, compatible with the load.
       * @param dofs Ordered map of local indices to prescribed coefficients.
       * @param offset Index shift applied to every key in @p dofs.
       * @return Reference to this system, with compressed operator storage.
       * @throws Alert::Exception If the operator is not square, the load has
       * the wrong size, or a shifted constraint index is out of bounds.
       */
      template <class DOFScalar>
      LinearSystem& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& stiffness = this->getOperator();
        auto& load = this->getVector();
        const size_t count = stiffness.cols();
        if (stiffness.rows() != stiffness.cols() || load.size() != stiffness.rows())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Elimination requires a square operator and a matching load."
            << Alert::Raise;
        }
        for (const auto& [index, value] : dofs)
        {
          if (offset >= count || index >= count - offset)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Shifted constraint index is outside the system." << Alert::Raise;
          }
        }
        stiffness.makeCompressed();
        if (dofs.empty())
          return *this;

        std::vector<unsigned char> constrained(count, 0);
        std::vector<SparseIndex> diagonals;
        diagonals.reserve(dofs.size());
        std::vector<SparseIndex> missing;
        auto* values = stiffness.valuePtr();
        auto* rows = stiffness.innerIndexPtr();
        auto* columns = stiffness.outerIndexPtr();
        auto* const rhs = load.data();
        for (const auto& [index, value] : dofs)
        {
          const SparseIndex column = index + offset;
          const DOFScalar prescribed = value;
          constrained[column] = 1;
          Optional<SparseIndex> diagonal;
          for (SparseIndex entry = columns[column]; entry < columns[column + 1]; ++entry)
          {
            rhs[rows[entry]] -= values[entry] * prescribed;
            if (rows[entry] == column)
              diagonal = entry;
          }
          if (diagonal)
            diagonals.push_back(*diagonal);
          else
            missing.push_back(column);
        }
        if (!missing.empty())
        {
          stiffness.resizeNonZeros(stiffness.nonZeros() + missing.size());
          values = stiffness.valuePtr();
          rows = stiffness.innerIndexPtr();
          columns = stiffness.outerIndexPtr();
          SparseIndex remaining = missing.size();
          for (SparseIndex column = count; column-- > 0;)
          {
            const SparseIndex start = columns[column];
            const SparseIndex end = columns[column + 1];
            SparseIndex source = end;
            SparseIndex destination = end + remaining;
            columns[column + 1] = destination;
            if (remaining && missing[remaining - 1] == column)
            {
              --remaining;
              while (source > start && rows[source - 1] > column)
              {
                --source;
                --destination;
                rows[destination] = rows[source];
                values[destination] = values[source];
              }
              --destination;
              rows[destination] = column;
              values[destination] = MatrixScalar(0);
            }
            while (source > start)
            {
              --source;
              --destination;
              rows[destination] = rows[source];
              values[destination] = values[source];
            }
          }
          diagonals.clear();
          for (const auto& [index, value] : dofs)
          {
            const SparseIndex column = index + offset;
            const auto* diagonal = std::lower_bound(
              rows + columns[column], rows + columns[column + 1], column);
            diagonals.push_back(diagonal - rows);
          }
        }

        constexpr size_t BlockSize = 16;
        constexpr size_t LocalitySamples = 16;
        constexpr size_t MinimumAverageColumnEntries = 8;
        const SparseIndex first = dofs.begin()->first + offset;
        const SparseIndex last = std::prev(dofs.end())->first + offset;
        bool prune = false;
        if (size_t(last - first + 1) == dofs.size() &&
          size_t(stiffness.nonZeros()) / count >= MinimumAverageColumnEntries)
        {
          size_t usable = 0;
          size_t rejected = 0;
          for (size_t sample = 0; sample < LocalitySamples; ++sample)
          {
            const SparseIndex position =
              size_t(stiffness.nonZeros()) * (2 * sample + 1) / (2 * LocalitySamples);
            const auto column =
              std::upper_bound(columns, columns + count + 1, position) - columns - 1;
            if (constrained[column])
              continue;
            ++usable;
            rejected +=
              rows[columns[column + 1] - 1] < first || rows[columns[column]] > last;
          }
          prune = usable && 2 * rejected >= usable;
        }

        SparseIndex begin = 0;
        auto constraint = dofs.begin();
        while (begin < stiffness.nonZeros())
        {
          const SparseIndex end = constraint == dofs.end()
            ? stiffness.nonZeros()
            : columns[constraint->first + offset];
          if (prune)
          {
            const auto firstColumn =
              std::upper_bound(columns, columns + count + 1, begin) - columns - 1;
            const auto lastColumn =
              constraint == dofs.end() ? count : constraint->first + offset;
            for (size_t column = firstColumn; column < lastColumn; ++column)
            {
              const SparseIndex start = columns[column];
              const SparseIndex stop = columns[column + 1];
              if (start == stop || rows[stop - 1] < first || rows[start] > last)
                continue;
              auto* lower = std::lower_bound(rows + start, rows + stop, first);
              const auto* upper = std::upper_bound(lower, rows + stop, last);
              std::fill(
                values + (lower - rows), values + (upper - rows), MatrixScalar(0));
            }
          }
          else
          {
            while (end - begin >= SparseIndex(BlockSize))
            {
              unsigned mask = 0;
              for (size_t bit = 0; bit < BlockSize; ++bit)
                mask |= unsigned(constrained[rows[begin + bit]]) << bit;
              while (mask)
              {
                values[begin + std::countr_zero(mask)] = MatrixScalar(0);
                mask &= mask - 1;
              }
              begin += BlockSize;
            }
            for (; begin < end; ++begin)
            {
              if constexpr (std::is_same_v<MatrixScalar, Real>)
              {
                const std::uint64_t mask = std::uint64_t(constrained[rows[begin]]) - 1;
                values[begin] =
                  std::bit_cast<Real>(std::bit_cast<std::uint64_t>(values[begin]) & mask);
              }
              else if (constrained[rows[begin]])
                values[begin] = MatrixScalar(0);
            }
          }
          if (constraint == dofs.end())
            break;
          SparseIndex column = constraint->first + offset;
          ++constraint;
          while (
            constraint != dofs.end() && constraint->first + offset == size_t(column + 1))
          {
            ++column;
            ++constraint;
          }
          begin = columns[column + 1];
          if (begin > end)
            std::fill(values + end, values + begin, MatrixScalar(0));
        }
        for (SparseIndex diagonal : diagonals)
          values[diagonal] = MatrixScalar(1);
        for (const auto& [index, value] : dofs)
          rhs[index + offset] = value;
        return *this;
      }

      /**
       * @brief Merges dependent degrees of freedom in the sparse linear system.
       *
       * Eliminates dependent DOFs by expressing them as linear combinations
       * of independent DOFs, modifying both the operator and RHS vector.
       * Column substitution uses the coefficients directly; restriction of
       * rows and the load uses their conjugates, giving the transformation
       * @f$T^*AT@f$ and @f$T^*b@f$ before restoring constraint rows.
       *
       * @tparam DOFScalar Scalar type for DOF values
       * @param[in] dofs Map from DOF index to (dependent indices, coefficients) pairs
       * @param[in] offset Index offset applied to DOF indices (default: 0)
       * @return Reference to @f$ *this @f$
       */
      template <class DOFScalar>
      LinearSystem& merge(
        const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs,
        size_t offset = 0)
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
          const Index localParent = q.back();
          const Index parent = localParent + static_cast<Index>(offset);
          assert(stiffness.rows() >= 0);
          assert(stiffness.cols() >= 0);
          assert(parent < static_cast<size_t>(stiffness.rows()));
          assert(parent < static_cast<size_t>(stiffness.cols()));
          q.pop_back();

          auto find = dofs.find(localParent);
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
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            stiffness.col(child) += coeff * stiffness.col(parent);
          }

          // Assumes CCS format
          for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(
                 stiffness, parent);
               it; ++it)
            it.valueRef() = 0;

          // Eliminate the parent row, adding it to the child rows
          IndexMap<MatrixScalar> parentLookup;
          std::vector<IndexMap<MatrixScalar>> childrenLookup(children.size());
          for (size_t col = 0; col < static_cast<size_t>(stiffness.cols()); col++)
          {
            Boolean parentFound = false;
            size_t childrenFound = 0;
            for (typename Math::SparseMatrix<MatrixScalar>::InnerIterator it(
                   stiffness, col);
                 it; ++it)
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
                    const Index child = children.coeff(i) + static_cast<Index>(offset);
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
              childrenLookup[i][col] += Math::conj(coeff) * value;
            }
          }

          for (size_t i = 0; i < count; i++)
          {
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            for (const auto& [col, value] : childrenLookup[i])
              stiffness.coeffRef(child, col) = value;
          }

          // Eliminate the parent entry, adding it to the child entries
          for (size_t i = 0; i < count; i++)
          {
            const MatrixScalar coeff = coeffs.coeff(i);
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            mass.coeffRef(child) += Math::conj(coeff) * mass.coeff(parent);
          }
          mass.coeffRef(parent) = 0;
        }

        for (const auto& [localParent, node] : dofs)
        {
          const Index parent = localParent + static_cast<Index>(offset);
          stiffness.coeffRef(parent, parent) = 1.0;
          const auto& [children, coeffs] = node;
          assert(children.size() >= 0);
          for (size_t i = 0; i < static_cast<size_t>(children.size()); i++)
          {
            const MatrixScalar coeff = coeffs.coeff(i);
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            stiffness.coeffRef(parent, child) = -coeff;
          }
        }
        return *this;
      }

      /**
       * @brief Returns the system operator matrix @f$ A @f$.
       * @returns The system operator matrix @f$ A @f$.
       */
      constexpr MatrixType& getOperator()
      {
        return m_operator;
      }

      /**
       * @brief Returns the system operator matrix @f$ A @f$.
       * @returns The system operator matrix @f$ A @f$.
       */
      constexpr const MatrixType& getOperator() const
      {
        return m_operator;
      }

      /**
       * @brief Returns the right-hand side vector @f$ b @f$.
       * @returns The right-hand side vector @f$ b @f$.
       */
      constexpr VectorType& getVector()
      {
        return m_vector;
      }

      /**
       * @brief Returns the right-hand side vector @f$ b @f$.
       * @returns The right-hand side vector @f$ b @f$.
       */
      constexpr const VectorType& getVector() const
      {
        return m_vector;
      }

      /**
       * @brief Returns the solution vector @f$ x @f$.
       * @returns The solution vector @f$ x @f$.
       */
      constexpr VectorType& getSolution()
      {
        return m_solution;
      }

      /**
       * @brief Returns the solution vector @f$ x @f$.
       * @returns The solution vector @f$ x @f$.
       */
      constexpr const VectorType& getSolution() const
      {
        return m_solution;
      }

    private:
      MatrixType m_operator; ///< The operator of the linear system.
      VectorType m_vector;   ///< The vector of the linear system.
      VectorType m_solution; ///< The solution vector of the linear system.
  };

  /**
   * @brief Linear system specialization for dense matrices.
   *
   * This specialization handles linear systems @f$ Ax = b @f$ where @f$ A @f$ is a 
   * dense matrix. This is typically used for smaller systems where the dense 
   * representation is more efficient.
   *
   * @tparam MatrixScalar Scalar type for matrix elements
   * @tparam VectorScalar Scalar type for vector elements
   */
  template <class MatrixScalar, class VectorScalar>
  class LinearSystem<Math::Matrix<MatrixScalar>, Math::Vector<VectorScalar>>
    : public LinearSystemBase<Math::Matrix<MatrixScalar>, Math::Vector<VectorScalar>,
        LinearSystem<Math::Matrix<MatrixScalar>, Math::Vector<VectorScalar>>>
  {
    public:
      /// Dense matrix type for the linear operator
      using MatrixType = Math::Matrix<MatrixScalar>;

      /// Vector type for solution and right-hand side
      using VectorType = Math::Vector<VectorScalar>;

      /// Parent class type
      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixType, VectorType>>;

      /// Default constructor
      constexpr LinearSystem() = default;

      /**
       * Copy constructor
       * @param other Object to copy from.
       */
      constexpr LinearSystem(const LinearSystem& other)
        : Parent(other),
          m_operator(other.m_operator),
          m_vector(other.m_vector),
          m_solution(other.m_solution)
      {}

      /**
       * Move constructor
       * @param other Object to move from.
       */
      constexpr LinearSystem(LinearSystem&& other) noexcept
        : Parent(std::move(other)),
          m_operator(std::move(other.m_operator)),
          m_vector(std::move(other.m_vector)),
          m_solution(std::move(other.m_solution))
      {}

      /// Virtual destructor
      virtual ~LinearSystem() = default;

      /**
       * @brief Copy assignment operator.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      constexpr LinearSystem& operator=(const LinearSystem& other)
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

      /**
       * @brief Move assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      constexpr LinearSystem& operator=(LinearSystem&& other) noexcept
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

      /**
       * @brief Eliminates degrees of freedom from the dense linear system.
       *
       * For each DOF in @p dofs, moves its contribution from LHS to RHS,
       * then zeros the corresponding row and column in the operator and
       * sets the diagonal to 1, effectively imposing the constraint.
       *
       * @tparam DOFScalar Scalar type for DOF values
       * @param[in] dofs Map from DOF index to prescribed value
       * @param[in] offset Index offset applied to DOF indices (default: 0)
       * @return Reference to @f$ *this @f$
       */
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

      /**
       * @brief Merges dependent degrees of freedom in the dense linear system.
       *
       * Eliminates dependent DOFs by expressing them as linear combinations
       * of independent DOFs, modifying both the operator and RHS vector.
       * Column substitution uses the coefficients directly; restriction of
       * rows and the load uses their conjugates, giving the transformation
       * @f$T^*AT@f$ and @f$T^*b@f$ before restoring constraint rows.
       *
       * @tparam DOFScalar Scalar type for DOF values
       * @param[in] dofs Map from DOF index to (dependent indices, coefficients) pairs
       * @param[in] offset Index offset applied to DOF indices (default: 0)
       * @return Reference to @f$ *this @f$
       */
      template <class DOFScalar>
      LinearSystem& merge(
        const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs,
        size_t offset = 0)
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
          const Index localParent = q.back();
          const Index parent = localParent + static_cast<Index>(offset);
          assert(stiffness.rows() >= 0);
          assert(stiffness.cols() >= 0);
          assert(parent < static_cast<size_t>(stiffness.rows()));
          assert(parent < static_cast<size_t>(stiffness.cols()));
          q.pop_back();

          auto find = dofs.find(localParent);
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
            const Index child = children.coeff(i) + static_cast<Index>(offset);
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
            for (typename Math::Matrix<MatrixScalar>::InnerIterator it(stiffness, col);
                 it; ++it)
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
                    const Index child = children.coeff(i) + static_cast<Index>(offset);
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
              childrenLookup[i][col] += Math::conj(coeff) * value;
            }
          }

          for (size_t i = 0; i < count; i++)
          {
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            for (const auto& [col, value] : childrenLookup[i])
              stiffness.coeffRef(child, col) = value;
          }

          // Eliminate the parent entry, adding it to the child entries
          for (size_t i = 0; i < count; i++)
          {
            const DOFScalar coeff = coeffs.coeff(i);
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            mass.coeffRef(child) += Math::conj(coeff) * mass.coeff(parent);
          }
          mass.coeffRef(parent) = 0;
        }

        for (const auto& [localParent, node] : dofs)
        {
          const Index parent = localParent + static_cast<Index>(offset);
          stiffness.coeffRef(parent, parent) = 1.0;
          const auto& [children, coeffs] = node;
          assert(children.size() >= 0);
          for (size_t i = 0; i < static_cast<size_t>(children.size()); i++)
          {
            const DOFScalar coeff = coeffs.coeff(i);
            const Index child = children.coeff(i) + static_cast<Index>(offset);
            stiffness.coeffRef(parent, child) = -coeff;
          }
        }
        return *this;
      }

      /**
       * @brief Returns the system operator matrix @f$ A @f$.
       * @returns The system operator matrix @f$ A @f$.
       */
      constexpr MatrixType& getOperator()
      {
        return m_operator;
      }

      /**
       * @brief Returns the system operator matrix @f$ A @f$.
       * @returns The system operator matrix @f$ A @f$.
       */
      constexpr const MatrixType& getOperator() const
      {
        return m_operator;
      }

      /**
       * @brief Returns the right-hand side vector @f$ b @f$.
       * @returns The right-hand side vector @f$ b @f$.
       */
      constexpr VectorType& getVector()
      {
        return m_vector;
      }

      /**
       * @brief Returns the right-hand side vector @f$ b @f$.
       * @returns The right-hand side vector @f$ b @f$.
       */
      constexpr const VectorType& getVector() const
      {
        return m_vector;
      }

      /**
       * @brief Returns the solution vector @f$ x @f$.
       * @returns The solution vector @f$ x @f$.
       */
      constexpr VectorType& getSolution()
      {
        return m_solution;
      }

      /**
       * @brief Returns the solution vector @f$ x @f$.
       * @returns The solution vector @f$ x @f$.
       */
      constexpr const VectorType& getSolution() const
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
    : std::integral_constant<std::size_t, 3>
  {};

  /// @brief Structured-binding element 0 of a LinearSystem: the operator @f$ A @f$.
  template <class Matrix, class Vector>
  struct tuple_element<0, Rodin::Math::LinearSystem<Matrix, Vector>>
  {
    /// @brief Type of element 0 (the operator matrix).
    using type = Matrix;
  };

  /// @brief Structured-binding element 1 of a LinearSystem: the solution @f$ x @f$.
  template <class Matrix, class Vector>
  struct tuple_element<1, Rodin::Math::LinearSystem<Matrix, Vector>>
  {
    /// @brief Type of element 1 (the solution vector).
    using type = Vector;
  };

  /// @brief Structured-binding element 2 of a LinearSystem: the right-hand side @f$ b @f$.
  template <class Matrix, class Vector>
  struct tuple_element<2, Rodin::Math::LinearSystem<Matrix, Vector>>
  {
    /// @brief Type of element 2 (the right-hand side vector).
    using type = Vector;
  };
}

namespace Rodin::Math
{
  /**
   * @brief Structured-binding accessor: element @p I of @p ls
   * (0 = operator @f$ A @f$, 1 = solution @f$ x @f$, 2 = right-hand side @f$ b @f$).
   * @param ls Linear system to solve.
   * @returns The .
   */
  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(LinearSystem<Matrix, Vector>& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0)
      return (ls.getOperator());
    else if constexpr (I == 1)
      return (ls.getSolution());
    else
      return (ls.getVector());
  }

  /**
   * @brief Structured-binding accessor (const overload); see get(LinearSystem&).
   * @param ls Linear system to solve.
   * @returns The .
   */
  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(const LinearSystem<Matrix, Vector>& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0)
      return (ls.getOperator());
    else if constexpr (I == 1)
      return (ls.getSolution());
    else
      return (ls.getVector());
  }

  /**
   * @brief Structured-binding accessor (rvalue overload); moves out element @p I.
   * @param ls Linear system to solve.
   * @returns The .
   */
  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(LinearSystem<Matrix, Vector>&& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0)
      return std::move(ls.getOperator());
    else if constexpr (I == 1)
      return std::move(ls.getSolution());
    else
      return std::move(ls.getVector());
  }

  /**
   * @brief Structured-binding accessor (const rvalue overload); see get(LinearSystem&&).
   * @param ls Linear system to solve.
   * @returns The .
   */
  template <std::size_t I, class Matrix, class Vector>
  decltype(auto) get(const LinearSystem<Matrix, Vector>&& ls)
  {
    static_assert(I < 3);
    if constexpr (I == 0)
      return std::move(ls.getOperator());
    else if constexpr (I == 1)
      return std::move(ls.getSolution());
    else
      return std::move(ls.getVector());
  }
}

#endif
