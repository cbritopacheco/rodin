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
 * shared with ParU through a non-owning CHOLMOD view. Rodin's default Eigen
 * sparse matrices use 32-bit indices whereas ParU requires 64-bit indices, so
 * only the two CSC index arrays are converted.
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
   * ParU factorizes a general matrix as @f$ PAQ = LU @f$. The symbolic and
   * numeric factorizations are rebuilt for every call to solve, matching the
   * assembly contract in which the matrix may change between solves.
   *
   * Architecture:
   * 1. Compress the Eigen CSC matrix and expose its values through
   *    Eigen::viewAsCholmod without copying them.
   * 2. Convert the CSC column pointers and row indices to the 64-bit layout
   *    required by ParU.
   * 3. Analyze, factorize, solve, and release all ParU objects with RAII.
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

      ParU(const ParU& other) = default;
      ParU(ParU&& other) noexcept = default;
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
        m_ordering = ordering;
        return *this;
      }

      /** @brief Returns the configured fill-reducing ordering. */
      Ordering getOrdering() const noexcept
      {
        return m_ordering;
      }

      /**
       * @brief Solves the linear system with ParU.
       * @param[in,out] axb System whose solution vector receives @f$ A^{-1}b @f$.
       */
      void solve(LinearSystemType& axb) override
      {
        auto& matrix = axb.getOperator();
        if (matrix.rows() != matrix.cols())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU requires a square matrix."
            << Alert::Raise;
        }
        if (axb.getVector().size() != matrix.rows())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix."
            << Alert::Raise;
        }

        matrix.makeCompressed();
        cholmod_sparse view = Eigen::viewAsCholmod(matrix);

        Array<std::int64_t> columnPointers(matrix.outerSize() + 1);
        Array<std::int64_t> rowIndices(matrix.nonZeros());
        for (Eigen::Index i = 0; i <= matrix.outerSize(); ++i)
          columnPointers(i) = static_cast<std::int64_t>(matrix.outerIndexPtr()[i]);
        for (Eigen::Index i = 0; i < matrix.nonZeros(); ++i)
          rowIndices(i) = static_cast<std::int64_t>(matrix.innerIndexPtr()[i]);

        view.p = columnPointers.data();
        view.i = rowIndices.data();
        view.itype = CHOLMOD_LONG;

        if (view.xtype != CHOLMOD_REAL || view.dtype != CHOLMOD_DOUBLE)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU requires a real double-precision CHOLMOD view, but got "
            << "xtype " << view.xtype << " and dtype " << view.dtype << "."
            << Alert::Raise;
        }

        Resources resources;
        check(ParU_InitControl(&resources.control), "control initialization");
        check(
          ParU_Set(
            PARU_CONTROL_MAX_THREADS,
            static_cast<std::int64_t>(m_maxThreads),
            resources.control),
          "thread-count configuration");
        if (m_ordering != Ordering::Default)
        {
          check(
            ParU_Set(
              PARU_CONTROL_ORDERING,
              static_cast<std::int64_t>(m_ordering),
              resources.control),
            "ordering configuration");
        }
        const auto analyzeInfo =
          ParU_Analyze(&view, &resources.symbolic, resources.control);
        if (analyzeInfo != PARU_SUCCESS)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ParU symbolic analysis failed with status "
            << static_cast<Integer>(analyzeInfo) << " for a "
            << view.nrow << " x " << view.ncol << " matrix (xtype "
            << view.xtype << ", dtype " << view.dtype << ")."
            << Alert::Raise;
        }
        check(
          ParU_Factorize(
            &view, resources.symbolic, &resources.numeric, resources.control),
          "numeric factorization");

        axb.getSolution().resize(axb.getVector().size());
        check(
          ParU_Solve(
            resources.symbolic,
            resources.numeric,
            axb.getVector().data(),
            axb.getSolution().data(),
            resources.control),
          "solve");
      }

      ParU* copy() const noexcept override
      {
        return new ParU(*this);
      }

    private:
      struct Resources
      {
        ParU_Control control = nullptr;
        ParU_Symbolic symbolic = nullptr;
        ParU_Numeric numeric = nullptr;

        ~Resources()
        {
          ParU_FreeNumeric(&numeric, control);
          ParU_FreeSymbolic(&symbolic, control);
          ParU_FreeControl(&control);
        }
      };

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
  };
}

#endif
#endif
