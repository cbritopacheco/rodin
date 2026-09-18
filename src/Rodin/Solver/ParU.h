/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ParU.h
 * @brief SuiteSparse:ParU parallel multifrontal sparse LU solver wrapper.
 *
 * ParU computes a sparse LU factorization of a real square matrix and uses
 * OpenMP tasks and parallel BLAS in its numerical phase. The matrix values are
 * shared with ParU through a non-owning CHOLMOD view. Rodin sparse matrices
 * use 32-bit indices whereas ParU requires 64-bit indices, so only the two
 * CSC index arrays are converted.
 *
 * @note This solver requires SuiteSparse:ParU and RODIN_USE_PARU.
 */
#ifndef RODIN_SOLVER_PARU_H
#define RODIN_SOLVER_PARU_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_PARU

#include <cstring>
#include <cstdint>
#include <limits>

#include <Eigen/CholmodSupport>
#include <ParU.h>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Array.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::ParU<LinearSystem>>
  {
      /// @brief Linear system type.
      using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup ParUSpecializations ParU Template Specializations
   * @brief Template specializations of the ParU class.
   * @see @ref ParU
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for ParU.
   */
  template <class LinearSystem, class TrialFunction, class TestFunction>
  ParU(Variational::Problem<LinearSystem, TrialFunction, TestFunction>&)
    -> ParU<LinearSystem>;

  /**
   * @ingroup ParUSpecializations
   * @brief Parallel multifrontal LU solver for real sparse square systems.
   *
   * ParU factorizes a general matrix as @f$ PAQ = LU @f$. Symbolic analysis is
   * reused while the matrix sparsity pattern and ordering remain unchanged.
   * The numeric factorization is retained after @ref factorize so multiple
   * right-hand sides can be solved without refactorization. The ordinary
   * @ref solve entry point compares the current matrix with the factorized
   * numeric state and refactorizes only when it has changed. Since Eigen does
   * not expose a matrix mutation counter, ParU retains an exact snapshot of
   * the factorized values for this comparison.
   *
   * Architecture:
   * 1. Compress the Eigen CSC matrix and expose its values through
   *    Eigen::viewAsCholmod without copying them.
   * 2. Convert the CSC column pointers and row indices to the signed 64-bit
   *    layout required by ParU.
   * 3. Reuse symbolic analysis when the converted CSC pattern is unchanged.
   * 4. Retain the numeric factorization for explicit repeated solves.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref ParU "ParU<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>>" | Real double-precision sparse systems. |
   */
  template <>
  class ParU<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>>
  {
    public:
      using ScalarType = Real;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::SparseMatrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /** @brief Fill-reducing ordering used during symbolic analysis. */
      enum class Ordering : std::int64_t
      {
        Default = -1,
        AMD = PARU_ORDERING_AMD,
        METIS = PARU_ORDERING_METIS,
        METISGuard = PARU_ORDERING_METIS_GUARD,
        CHOLMOD = PARU_ORDERING_CHOLMOD,
        Best = PARU_ORDERING_BEST,
        Natural = PARU_ORDERING_NONE
      };

      /** @brief Constructs a ParU solver for the given problem. */
      ParU(ProblemBaseType& pb)
        : Parent(pb)
      {}

      ParU(const ParU& other)
        : Parent(other),
          m_maxThreads(other.m_maxThreads),
          m_ordering(other.m_ordering)
      {}

      ParU(ParU&& other) noexcept
        : Parent(std::move(other)),
          m_maxThreads(other.m_maxThreads),
          m_ordering(other.m_ordering)
      {}

      ~ParU() override = default;

      /**
       * @brief Sets ParU's maximum thread count.
       * @param count Maximum threads, or zero to use ParU's OpenMP default.
       * @returns Reference to this solver.
       */
      ParU& setMaxThreads(Index count)
      {
        if (count > static_cast<Index>(std::numeric_limits<std::int64_t>::max()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The maximum thread count exceeds ParU's integer range."
            << Alert::Raise;
        }
        m_maxThreads = count;
        return *this;
      }

      /** @brief Returns ParU's configured maximum thread count. */
      Index getMaxThreads() const noexcept
      {
        return m_maxThreads;
      }

      /**
       * @brief Sets the fill-reducing ordering.
       * @param ordering Ordering strategy, or Ordering::Default for ParU's default.
       * @returns Reference to this solver.
       */
      ParU& setOrdering(Ordering ordering) noexcept
      {
        if (m_ordering != ordering)
        {
          m_ordering = ordering;
          m_resources.clearSymbolic();
          m_factorizedValues.resize(0);
        }
        return *this;
      }

      /** @brief Returns the configured fill-reducing ordering. */
      Ordering getOrdering() const noexcept
      {
        return m_ordering;
      }

      /**
       * @brief Analyzes and numerically factorizes the system matrix.
       *
       * Reuses symbolic analysis when the sparsity pattern and ordering match
       * the previous factorization. Any previous numeric factorization is
       * replaced.
       */
      void factorize(LinearSystemType& axb)
      {
        auto& matrix = axb.getOperator();
        if (matrix.rows() != matrix.cols())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU requires a square matrix."
            << Alert::Raise;
        }

        matrix.makeCompressed();
        cholmod_sparse view = Eigen::viewAsCholmod(matrix);

        m_resources.initialize(*this);
        const bool samePattern = hasSamePattern(matrix);
        if (!samePattern)
        {
          m_resources.clearSymbolic();
          m_columnPointers.resize(matrix.outerSize() + 1);
          m_rowIndices.resize(matrix.nonZeros());
          for (Eigen::Index i = 0; i <= matrix.outerSize(); ++i)
          {
            m_columnPointers(i) =
              static_cast<SuiteSparse_long>(matrix.outerIndexPtr()[i]);
          }
          for (Eigen::Index i = 0; i < matrix.nonZeros(); ++i)
          {
            m_rowIndices(i) =
              static_cast<SuiteSparse_long>(matrix.innerIndexPtr()[i]);
          }
        }

        view.p = m_columnPointers.data();
        view.i = m_rowIndices.data();
        view.itype = CHOLMOD_LONG;

        if (view.xtype != CHOLMOD_REAL || view.dtype != CHOLMOD_DOUBLE)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU requires a real double-precision CHOLMOD view, but got "
            << "xtype " << view.xtype << " and dtype " << view.dtype << "."
            << Alert::Raise;
        }

        configureControl();
        if (!m_resources.symbolic)
        {
          const auto analyzeInfo =
            ParU_Analyze(
              &view, &m_resources.symbolic, m_resources.control);
          if (analyzeInfo != PARU_SUCCESS)
          {
            m_resources.clearSymbolic();
            Alert::MemberFunctionException(*this, __func__)
              << "ParU symbolic analysis failed with status "
              << static_cast<Integer>(analyzeInfo) << " for a "
              << view.nrow << " x " << view.ncol << " matrix (xtype "
              << view.xtype << ", dtype " << view.dtype << ")."
              << Alert::Raise;
          }
        }

        m_resources.clearNumeric();
        m_factorizedValues.resize(0);
        const auto factorizeInfo =
          ParU_Factorize(
            &view,
            m_resources.symbolic,
            &m_resources.numeric,
            m_resources.control);
        if (factorizeInfo != PARU_SUCCESS)
        {
          m_resources.clearNumeric();
          check(factorizeInfo, "numeric factorization");
        }

        m_factorizedValues.resize(matrix.nonZeros());
        for (Eigen::Index i = 0; i < matrix.nonZeros(); ++i)
          m_factorizedValues(i) = matrix.valuePtr()[i];
      }

      /**
       * @brief Solves using the retained numeric factorization.
       *
       * The caller must call @ref factorize again after changing the matrix.
       * Only the right-hand side may change between calls.
       */
      void solveFactorized(LinearSystemType& axb)
      {
        if (!m_resources.numeric)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No ParU numeric factorization is available."
            << Alert::Raise;
        }
        if (axb.getVector().size() != axb.getOperator().rows())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix."
            << Alert::Raise;
        }

        configureControl();
        axb.getSolution().resize(axb.getVector().size());
        check(
          ParU_Solve(
            m_resources.symbolic,
            m_resources.numeric,
            axb.getVector().data(),
            axb.getSolution().data(),
            m_resources.control),
          "solve");
      }

      /** @brief Returns whether a reusable numeric factorization is available. */
      bool hasFactorization() const noexcept
      {
        return m_resources.numeric != nullptr;
      }

      /**
       * @brief Returns whether the retained factorization matches a system.
       *
       * Both the CSC sparsity pattern and the matrix values are compared.
       */
      bool hasFactorization(const LinearSystemType& axb) const noexcept
      {
        const auto& matrix = axb.getOperator();
        if (!m_resources.numeric || !matrix.isCompressed() ||
            !hasSamePattern(matrix) ||
            m_factorizedValues.size() != matrix.nonZeros())
          return false;

        if (matrix.nonZeros() == 0)
          return true;

        return std::memcmp(
          m_factorizedValues.data(),
          matrix.valuePtr(),
          static_cast<size_t>(matrix.nonZeros()) * sizeof(ScalarType)) == 0;
      }

      /** @brief Releases the retained numeric factorization. */
      void clearFactorization() noexcept
      {
        m_resources.clearNumeric();
        m_factorizedValues.resize(0);
      }

      /**
       * @brief Solves the linear system with ParU.
       * @param[in,out] axb System whose solution vector receives @f$ A^{-1}b @f$.
       */
      void solve(LinearSystemType& axb) override
      {
        if (!hasFactorization(axb))
          factorize(axb);
        solveFactorized(axb);
      }

      ParU* copy() const noexcept override
      {
        return new ParU(*this);
      }

    private:
      struct Resources
      {
        Resources() = default;
        Resources(const Resources&) = delete;
        Resources& operator=(const Resources&) = delete;

        ParU_Control control = nullptr;
        ParU_Symbolic symbolic = nullptr;
        ParU_Numeric numeric = nullptr;

        void initialize(const ParU& solver)
        {
          if (!control)
            solver.check(ParU_InitControl(&control), "control initialization");
        }

        void clearNumeric() noexcept
        {
          if (numeric)
            ParU_FreeNumeric(&numeric, control);
        }

        void clearSymbolic() noexcept
        {
          clearNumeric();
          if (symbolic)
            ParU_FreeSymbolic(&symbolic, control);
        }

        ~Resources()
        {
          clearSymbolic();
          if (control)
            ParU_FreeControl(&control);
        }
      };

      void configureControl()
      {
        check(
          ParU_Set(
            PARU_CONTROL_MAX_THREADS,
            static_cast<std::int64_t>(m_maxThreads),
            m_resources.control),
          "thread-count configuration");
        check(
          ParU_Set(
            PARU_CONTROL_ORDERING,
            m_ordering == Ordering::Default
              ? static_cast<std::int64_t>(PARU_DEFAULT_ORDERING)
              : static_cast<std::int64_t>(m_ordering),
            m_resources.control),
          "ordering configuration");
      }

      bool hasSamePattern(const OperatorType& matrix) const
      {
        if (m_columnPointers.size() != matrix.outerSize() + 1 ||
            m_rowIndices.size() != matrix.nonZeros())
          return false;

        for (Eigen::Index i = 0; i <= matrix.outerSize(); ++i)
        {
          if (m_columnPointers(i) != matrix.outerIndexPtr()[i])
            return false;
        }
        for (Eigen::Index i = 0; i < matrix.nonZeros(); ++i)
        {
          if (m_rowIndices(i) != matrix.innerIndexPtr()[i])
            return false;
        }
        return true;
      }

      void check(ParU_Info info, StringView operation) const
      {
        if (info != PARU_SUCCESS)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU " << operation << " failed with status "
            << static_cast<Integer>(info) << "."
            << Alert::Raise;
        }
      }

      Index m_maxThreads = 0;
      Ordering m_ordering = Ordering::Default;
      Array<SuiteSparse_long> m_columnPointers;
      Array<SuiteSparse_long> m_rowIndices;
      Array<ScalarType> m_factorizedValues;
      Resources m_resources;
  };
}

#endif
#endif
