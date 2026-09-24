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
#include "Info.h"
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
  template <class LinearSystem>
  ParU(Variational::ProblemBase<LinearSystem>&) -> ParU<LinearSystem>;

  /**
   * @ingroup ParUSpecializations
   * @brief Parallel multifrontal LU solver for real sparse square systems.
   *
   * ParU factorizes a general matrix as @f$ PAQ = LU @f$. Symbolic analysis is
   * reused across explicit @ref factorize calls, which assume that the matrix
   * sparsity pattern is unchanged.
   * The numeric factorization is retained after @ref factorize so multiple
   * right-hand sides can be solved without refactorization. The ordinary
   * @ref solve entry point builds a factorization only when none is retained.
   * Call @ref factorize explicitly after changing matrix values. After changing
   * the sparsity pattern, clear the symbolic factorization so the next
   * factorization also repeats symbolic analysis.
   *
   * Architecture:
   * 1. Compress the Eigen CSC matrix and expose its values through
   *    Eigen::viewAsCholmod without copying them.
   * 2. Convert the CSC column pointers and row indices to the signed 64-bit
   *    layout required by ParU.
   * 3. Reuse symbolic analysis until explicitly cleared.
   * 4. Retain the numeric factorization for explicit repeated solves.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref ParU "ParU<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>>" | Real double-precision sparse systems. |
   */
  template <>
  class ParU<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>> final
    : public LinearSolverBase<
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>>
  {
    public:
      using ScalarType = Real;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::SparseMatrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /** @brief Factorization stage, shared by the factorization solvers. */
      using Factorization = Rodin::Solver::Factorization;

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

      /** @brief Factorization stages retained between solves. */
      struct Resources
      {
          Resources() = default;
          Resources(const Resources&) = delete;
          Resources& operator=(const Resources&) = delete;

          /// @brief ParU control object, or nullptr before initialization.
          ParU_Control control = nullptr;
          /// @brief Symbolic analysis, or nullptr when none is retained.
          ParU_Symbolic symbolic = nullptr;
          /// @brief Numeric factorization, or nullptr when none is retained.
          ParU_Numeric numeric = nullptr;

          void initialize(const ParU& solver)
          {
            if (!control)
              solver.check(ParU_InitControl(&control), "control initialization");
          }

          /**
           * @brief Releases a factorization stage.
           *
           * Releasing Factorization::Symbolic also releases its dependent
           * numeric factorization.
           */
          void clear(Factorization factorization) noexcept
          {
            if (numeric)
              ParU_FreeNumeric(&numeric, control);
            if (factorization == Factorization::Symbolic && symbolic)
              ParU_FreeSymbolic(&symbolic, control);
          }

          ~Resources()
          {
            clear(Factorization::Symbolic);
            if (control)
              ParU_FreeControl(&control);
          }
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
            << "The maximum thread count exceeds ParU's integer range." << Alert::Raise;
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
          clear(Factorization::Symbolic);
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
       * Reuses existing symbolic analysis and replaces any previous numeric
       * factorization. The matrix sparsity pattern must remain unchanged until
       * the symbolic factorization is cleared.
       */
      void factorize(LinearSystemType& axb)
      {
        auto& matrix = axb.getOperator();
        if (matrix.rows() != matrix.cols())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU requires a square matrix." << Alert::Raise;
        }

        matrix.makeCompressed();
        cholmod_sparse view = Eigen::viewAsCholmod(matrix);

        m_resources.initialize(*this);
        if (!m_resources.symbolic)
        {
          m_columnPointers.resize(matrix.outerSize() + 1);
          m_rowIndices.resize(matrix.nonZeros());
          for (Eigen::Index i = 0; i <= matrix.outerSize(); ++i)
          {
            m_columnPointers(i) =
              static_cast<SuiteSparse_long>(matrix.outerIndexPtr()[i]);
          }
          for (Eigen::Index i = 0; i < matrix.nonZeros(); ++i)
          {
            m_rowIndices(i) = static_cast<SuiteSparse_long>(matrix.innerIndexPtr()[i]);
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

        configure();
        if (!m_resources.symbolic)
        {
          const auto analyzeInfo =
            ParU_Analyze(&view, &m_resources.symbolic, m_resources.control);
          if (analyzeInfo != PARU_SUCCESS)
          {
            m_resources.clear(Factorization::Symbolic);
            return record(analyzeInfo, {});
          }
        }
        m_info.factorization = Factorization::Symbolic;

        m_resources.clear(Factorization::Numeric);
        const auto factorizeInfo = ParU_Factorize(
          &view, m_resources.symbolic, &m_resources.numeric, m_resources.control);
        if (factorizeInfo != PARU_SUCCESS)
        {
          m_resources.clear(Factorization::Numeric);
          return record(factorizeInfo, Factorization::Symbolic);
        }
        record(PARU_SUCCESS, Factorization::Numeric);
      }

      /**
       * @brief Returns the outcome of the most recent operation.
       *
       * Updated by every factorization and every solve, so a caller reads it
       * after the call it wants to check.
       */
      const Info& getInfo() const noexcept
      {
        return m_info;
      }

      /**
       * @brief Checks whether the most recent operation succeeded.
       * @returns true if the solver succeeded, false otherwise.
       */
      Boolean success() const noexcept
      {
        return m_info.success;
      }

      /**
       * @brief Returns the retained ParU resources.
       *
       * A non-null @ref Resources::numeric denotes a reusable numeric
       * factorization, and a non-null @ref Resources::symbolic denotes reusable
       * symbolic analysis.
       */
      const Resources& getResources() const noexcept
      {
        return m_resources;
      }

      /**
       * @brief Releases a retained factorization stage.
       *
       * Clearing Factorization::Numeric preserves symbolic analysis. Clearing
       * Factorization::Symbolic also releases its dependent numeric
       * factorization and the converted CSC structure.
       */
      void clear(Factorization factorization) noexcept
      {
        m_resources.clear(factorization);
        if (factorization == Factorization::Symbolic)
        {
          m_columnPointers.resize(0);
          m_rowIndices.resize(0);
          m_info.factorization.reset();
        }
        else if (m_resources.symbolic)
        {
          m_info.factorization = Factorization::Symbolic;
        }
        else
        {
          m_info.factorization.reset();
        }
      }

      /**
       * @brief Solves the linear system with ParU.
       *
       * Reuses the retained numeric factorization and builds one only when
       * none is retained, so repeated right-hand sides are solved without
       * refactorization. The caller must call @ref factorize after changing
       * the matrix.
       *
       * @param[in,out] axb System whose solution vector receives @f$ A^{-1}b @f$.
       */
      void solve(LinearSystemType& axb) override
      {
        if (axb.getVector().size() != axb.getOperator().rows())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }

        if (!getResources().numeric)
        {
          factorize(axb);
          if (!success())
            return;
        }

        configure();
        axb.getSolution().resize(axb.getVector().size());
        record(ParU_Solve(m_resources.symbolic, m_resources.numeric,
                 axb.getVector().data(), axb.getSolution().data(), m_resources.control),
          Factorization::Numeric);
      }

      ParU* copy() const noexcept override
      {
        return new ParU(*this);
      }

    private:
      void configure()
      {
        check(ParU_Set(PARU_CONTROL_MAX_THREADS, static_cast<std::int64_t>(m_maxThreads),
                m_resources.control),
          "thread-count configuration");
        check(ParU_Set(PARU_CONTROL_ORDERING,
                m_ordering == Ordering::Default
                  ? static_cast<std::int64_t>(PARU_DEFAULT_ORDERING)
                  : static_cast<std::int64_t>(m_ordering),
                m_resources.control),
          "ordering configuration");
      }

      /**
       * @brief Records a ParU status and the stage that survived it.
       *
       * A failing analysis, factorization or solve is an outcome rather than a
       * defect in the caller, so it is reported through getInfo() instead of
       * raising.
       */
      void record(ParU_Info info, Optional<Factorization> factorization)
      {
        m_info.status = static_cast<Integer>(info);
        m_info.success = info == PARU_SUCCESS;
        m_info.factorization = factorization;
      }

      void check(ParU_Info info, StringView operation) const
      {
        if (info != PARU_SUCCESS)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU " << operation << " failed with status "
            << static_cast<Integer>(info) << "." << Alert::Raise;
        }
      }

      /// @brief Outcome of the most recent operation.
      Info m_info;

      Index m_maxThreads = 0;
      Ordering m_ordering = Ordering::Default;
      Array<SuiteSparse_long> m_columnPointers;
      Array<SuiteSparse_long> m_rowIndices;
      Resources m_resources;
  };
}

#endif
#endif
