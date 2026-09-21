/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MUMPS.h
 * @brief MUMPS multifrontal sparse direct solver wrapper.
 *
 * MUMPS factorizes a real square matrix with a multifrontal method, in
 * parallel over MPI processes and threaded BLAS. A symmetric matrix is
 * factorized as @f$ LDL^T @f$, which stores roughly half the entries of the
 * corresponding @f$ LU @f$ factorization.
 *
 * MUMPS reads the matrix in coordinate format with one-based indices, whereas
 * Rodin stores it as zero-based CSC. Only the two index arrays are converted,
 * once per sparsity pattern. The values of an unsymmetric matrix are shared
 * with MUMPS without a copy; a symmetric matrix is passed as its lower
 * triangle and therefore requires one value copy per factorization.
 *
 * @note This solver requires MUMPS and RODIN_USE_MUMPS.
 */
#ifndef RODIN_SOLVER_MUMPS_H
#define RODIN_SOLVER_MUMPS_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_MUMPS

#include <limits>

#include <dmumps_c.h>
#include <mpi.h>

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
  struct Traits<Solver::MUMPS<LinearSystem>>
  {
      /// @brief Linear system type.
      using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup MUMPSSpecializations MUMPS Template Specializations
   * @brief Template specializations of the MUMPS class.
   * @see @ref MUMPS
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for MUMPS.
   */
  template <class LinearSystem>
  MUMPS(Variational::ProblemBase<LinearSystem>&) -> MUMPS<LinearSystem>;

  /**
   * @ingroup MUMPSSpecializations
   * @brief Multifrontal direct solver for real sparse square systems.
   *
   * MUMPS factorizes a general matrix as @f$ A = LU @f$ and a symmetric matrix
   * as @f$ A = LDL^T @f$. Symbolic analysis is reused across explicit
   * @ref factorize calls, which assume that the matrix sparsity pattern is
   * unchanged. The numeric factorization is retained after @ref factorize so
   * multiple right-hand sides can be solved without refactorization. The
   * ordinary @ref solve entry point builds a factorization only when none is
   * retained. Call @ref factorize explicitly after changing matrix values.
   * After changing the sparsity pattern, clear the symbolic factorization so
   * the next factorization also repeats symbolic analysis.
   *
   * Symmetry is never inferred from the entries: a matrix is factorized as
   * symmetric only after @ref setSymmetric, which selects the triangle MUMPS
   * reads and the factorization it performs. Declaring a matrix symmetric
   * halves the stored factor, so it is worth doing whenever the assembled
   * operator is symmetric.
   *
   * Architecture:
   * 1. Convert the zero-based CSC structure to one-based coordinate arrays,
   *    keeping only the lower triangle for a symmetric matrix.
   * 2. Share the values of an unsymmetric matrix with MUMPS without a copy.
   * 3. Reuse symbolic analysis until explicitly cleared.
   * 4. Retain the numeric factorization for explicit repeated solves.
   *
   * MUMPS is parallel over MPI processes, and within each process over OpenMP
   * threads and threaded BLAS. This specialization solves on @c MPI_COMM_SELF,
   * so its parallelism and its memory are both shared: the whole factorization
   * is held in this process and computed by @ref setMaxThreads threads. Under
   * a distributed run each process solves its own system independently.
   *
   * MUMPS calls into MPI. When the host application has not initialized MPI,
   * the first factorization initializes it and leaves finalization to the
   * host, so that a solver instance never outlives the MPI session it uses.
   *
   * The system matrix must outlive the solves that follow a factorization,
   * because its values may be shared with MUMPS rather than copied.
   *
   * Info::status holds MUMPS's own @c INFOG(1): zero on success, a negative
   * value for an error, and a positive one for a warning.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref MUMPS "MUMPS<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>>" | Real double-precision sparse systems. |
   */
  template <>
  class MUMPS<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>> final
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

      /** @brief Matrix symmetry declared to MUMPS. */
      enum class Symmetry : int
      {
        /// @brief No symmetry; the whole matrix is factorized as @f$ LU @f$.
        Unsymmetric = 0,
        /// @brief Symmetric positive definite; factorized as @f$ LDL^T @f$
        /// without pivoting.
        PositiveDefinite = 1,
        /// @brief Symmetric, possibly indefinite; factorized as @f$ LDL^T @f$.
        General = 2
      };

      /** @brief Fill-reducing ordering used during symbolic analysis. */
      enum class Ordering : int
      {
        AMD = 0,
        AMF = 2,
        Scotch = 3,
        PORD = 4,
        METIS = 5,
        QAMD = 6,
        /// @brief Let MUMPS choose among the orderings it was built with.
        Automatic = 7
      };

      /** @brief Factorization stages retained between solves. */
      struct Resources
      {
          Resources() = default;
          Resources(const Resources&) = delete;
          Resources& operator=(const Resources&) = delete;

          /// @brief MUMPS instance, valid while @ref initialized is true.
          DMUMPS_STRUC_C instance{};
          /// @brief Whether the MUMPS instance has been initialized.
          bool initialized = false;
          /// @brief Whether reusable symbolic analysis is retained.
          bool symbolic = false;
          /// @brief Whether a reusable numeric factorization is retained.
          bool numeric = false;

          void initialize(const MUMPS& solver, Symmetry symmetry)
          {
            if (initialized)
              return;
            // This specialization is nailed to the process-local linear system
            // type, which only Local-context assembly produces, so the solve is
            // private to this process. MPI_COMM_WORLD would instead make it
            // collective over ranks that each hold a different matrix. This
            // mirrors the PETSc backend, which creates its objects on
            // PETSC_COMM_SELF for a Context::Local mesh.
            instance.comm_fortran = static_cast<MUMPS_INT>(MPI_Comm_c2f(MPI_COMM_SELF));
            instance.par = 1;
            instance.sym = static_cast<MUMPS_INT>(symmetry);
            instance.job = JobInitialize;
            dmumps_c(&instance);
            initialized = true;
            solver.check("initialization");
          }

          /**
           * @brief Releases retained factorization resources.
           *
           * MUMPS resources are released coarsely through destruction of the
           * current instance. Consequently, clearing either retained stage
           * destroys the MUMPS instance and drops any dependent retained
           * state. Releasing Factorization::Symbolic also releases its
           * dependent numeric factorization.
           */
          void clear(Factorization factorization) noexcept
          {
            if (!initialized)
            {
              numeric = false;
              if (factorization == Factorization::Symbolic)
                symbolic = false;
              return;
            }

            if (factorization == Factorization::Numeric)
            {
              if (!numeric)
                return;
            }
            else
            {
              if (!symbolic && !numeric)
                return;
            }

            destroy();
          }

          void destroy() noexcept
          {
            numeric = false;
            symbolic = false;
            if (!initialized)
              return;
            instance.job = JobDestroy;
            dmumps_c(&instance);
            initialized = false;
          }

          ~Resources()
          {
            destroy();
          }
      };

      /** @brief Constructs a MUMPS solver for the given problem. */
      MUMPS(ProblemBaseType& pb)
        : Parent(pb)
      {}

      MUMPS(const MUMPS& other)
        : Parent(other),
          m_maxThreads(other.m_maxThreads),
          m_symmetry(other.m_symmetry),
          m_ordering(other.m_ordering),
          m_workspacePercentage(other.m_workspacePercentage)
      {}

      MUMPS(MUMPS&& other) noexcept
        : Parent(std::move(other)),
          m_maxThreads(other.m_maxThreads),
          m_symmetry(other.m_symmetry),
          m_ordering(other.m_ordering),
          m_workspacePercentage(other.m_workspacePercentage)
      {}

      ~MUMPS() override = default;

      /**
       * @brief Declares the symmetry of the system matrix.
       *
       * The caller guarantees that the matrix is symmetric for
       * Symmetry::PositiveDefinite and Symmetry::General: only its lower
       * triangle is read. Symmetry selects the factorization itself and
       * therefore releases any retained factorization.
       *
       * @param symmetry Symmetry of the matrices to be factorized.
       * @returns Reference to this solver.
       */
      MUMPS& setSymmetric(Symmetry symmetry) noexcept
      {
        if (m_symmetry != symmetry)
        {
          m_symmetry = symmetry;
          clear(Factorization::Symbolic);
          m_resources.destroy();
          m_info = Info{};
        }
        return *this;
      }

      /** @brief Returns the declared symmetry of the system matrix. */
      Symmetry getSymmetric() const noexcept
      {
        return m_symmetry;
      }

      /**
       * @brief Sets the number of OpenMP threads MUMPS may use.
       * @param count Maximum threads, or zero to use the OpenMP default.
       * @returns Reference to this solver.
       */
      MUMPS& setMaxThreads(Index count)
      {
        if (count > static_cast<Index>(std::numeric_limits<MUMPS_INT>::max()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The maximum thread count exceeds the MUMPS integer range."
            << Alert::Raise;
        }
        m_maxThreads = count;
        return *this;
      }

      /** @brief Returns the configured maximum thread count. */
      Index getMaxThreads() const noexcept
      {
        return m_maxThreads;
      }

      /**
       * @brief Sets the fill-reducing ordering.
       * @param ordering Ordering strategy, or Ordering::Automatic to let MUMPS
       * choose.
       * @returns Reference to this solver.
       */
      MUMPS& setOrdering(Ordering ordering) noexcept
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
       * @brief Sets the working space MUMPS may add beyond its estimate.
       * @param percentage Percentage increase, or zero for the MUMPS default.
       * @returns Reference to this solver.
       */
      MUMPS& setWorkspacePercentage(Integer percentage)
      {
        if (percentage < 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The working space percentage must not be negative." << Alert::Raise;
        }
        m_workspacePercentage = percentage;
        return *this;
      }

      /** @brief Returns the configured working space percentage. */
      Integer getWorkspacePercentage() const noexcept
      {
        return m_workspacePercentage;
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
            << "MUMPS requires a square matrix." << Alert::Raise;
        }
        if (matrix.rows() > std::numeric_limits<MUMPS_INT>::max())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The matrix has " << matrix.rows()
            << " rows, which exceeds the MUMPS integer range." << Alert::Raise;
        }

        matrix.makeCompressed();
        initializeSession();
        m_resources.initialize(*this, m_symmetry);

        if (!m_resources.symbolic)
          buildIndices(matrix);

        auto& instance = m_resources.instance;
        instance.n = static_cast<MUMPS_INT>(matrix.rows());
        instance.nnz = static_cast<MUMPS_INT8>(m_rowIndices.size());
        instance.irn = m_rowIndices.data();
        instance.jcn = m_columnIndices.data();
        if (m_symmetry == Symmetry::Unsymmetric)
        {
          // The values of an unsymmetric matrix are read in place.
          instance.a = const_cast<ScalarType*>(matrix.valuePtr());
        }
        else
        {
          copyLowerTriangle(matrix);
          instance.a = m_values.data();
        }

        configure();
        if (!m_resources.symbolic)
        {
          instance.job = JobAnalyze;
          dmumps_c(&instance);
          if (instance.infog[0] < 0)
          {
            m_resources.clear(Factorization::Symbolic);
            return record({});
          }
          m_resources.symbolic = true;
        }

        m_resources.clear(Factorization::Numeric);
        instance.job = JobFactorize;
        dmumps_c(&instance);
        if (instance.infog[0] < 0)
          return record(Factorization::Symbolic);
        m_resources.numeric = true;
        record(Factorization::Numeric);
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
       * @brief Returns the retained MUMPS resources.
       *
       * A true @ref Resources::numeric denotes a reusable numeric
       * factorization, and a true @ref Resources::symbolic denotes reusable
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
       * factorization and the converted coordinate structure.
       */
      void clear(Factorization factorization) noexcept
      {
        m_resources.clear(factorization);
        if (factorization == Factorization::Symbolic)
        {
          m_rowIndices.resize(0);
          m_columnIndices.resize(0);
          m_values.resize(0);
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
       * @brief Solves the linear system with MUMPS.
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

        // MUMPS solves in place, overwriting the right-hand side it is given.
        axb.getSolution() = axb.getVector();
        auto& instance = m_resources.instance;
        instance.nrhs = 1;
        instance.lrhs = static_cast<MUMPS_INT>(axb.getSolution().size());
        instance.rhs = axb.getSolution().data();
        instance.job = JobSolve;
        dmumps_c(&instance);
        record(Factorization::Numeric);
      }

      MUMPS* copy() const noexcept override
      {
        return new MUMPS(*this);
      }

    private:
      static constexpr MUMPS_INT JobInitialize = -1;
      static constexpr MUMPS_INT JobDestroy = -2;
      static constexpr MUMPS_INT JobAnalyze = 1;
      static constexpr MUMPS_INT JobFactorize = 2;
      static constexpr MUMPS_INT JobSolve = 3;

      /**
       * @brief Ensures that MPI is available to MUMPS.
       *
       * Finalization is left to the host application: a solver instance
       * destroyed after @c MPI_Finalize could not release its MUMPS instance.
       */
      void initializeSession() const
      {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (finalized)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MPI has been finalized, so MUMPS can no longer be used." << Alert::Raise;
        }

        int initialized = 0;
        MPI_Initialized(&initialized);
        if (initialized)
          return;

        int argc = 0;
        char** argv = nullptr;
        int provided = 0;
        if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MPI initialization failed." << Alert::Raise;
        }
      }

      /**
       * @brief Converts the CSC structure to one-based coordinate arrays.
       *
       * A symmetric matrix contributes only its lower triangle, which is the
       * triangle MUMPS reads.
       */
      void buildIndices(const OperatorType& matrix)
      {
        const auto* outer = matrix.outerIndexPtr();
        const auto* inner = matrix.innerIndexPtr();
        const Eigen::Index columns = matrix.outerSize();
        const bool symmetric = m_symmetry != Symmetry::Unsymmetric;

        Eigen::Index count = 0;
        for (Eigen::Index column = 0; column < columns; ++column)
        {
          for (Eigen::Index k = outer[column]; k < outer[column + 1]; ++k)
          {
            if (symmetric && inner[k] < column)
              continue;
            ++count;
          }
        }

        m_rowIndices.resize(count);
        m_columnIndices.resize(count);
        Eigen::Index entry = 0;
        for (Eigen::Index column = 0; column < columns; ++column)
        {
          for (Eigen::Index k = outer[column]; k < outer[column + 1]; ++k)
          {
            if (symmetric && inner[k] < column)
              continue;
            m_rowIndices(entry) = static_cast<MUMPS_INT>(inner[k]) + 1;
            m_columnIndices(entry) = static_cast<MUMPS_INT>(column) + 1;
            ++entry;
          }
        }
      }

      /**
       * @brief Copies the lower triangle values in the order of @ref buildIndices.
       */
      void copyLowerTriangle(const OperatorType& matrix)
      {
        const auto* outer = matrix.outerIndexPtr();
        const auto* inner = matrix.innerIndexPtr();
        const auto* values = matrix.valuePtr();
        const Eigen::Index columns = matrix.outerSize();

        m_values.resize(m_rowIndices.size());
        Eigen::Index entry = 0;
        for (Eigen::Index column = 0; column < columns; ++column)
        {
          for (Eigen::Index k = outer[column]; k < outer[column + 1]; ++k)
          {
            if (inner[k] < column)
              continue;
            m_values(entry) = values[k];
            ++entry;
          }
        }

        if (entry != m_values.size())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The matrix sparsity pattern changed since the symbolic analysis."
            << Alert::Raise;
        }
      }

      void configure()
      {
        auto& instance = m_resources.instance;
        // Silence the MUMPS error, warning, and diagnostic streams.
        instance.icntl[0] = -1;
        instance.icntl[1] = -1;
        instance.icntl[2] = -1;
        instance.icntl[3] = 0;
        instance.icntl[6] = static_cast<MUMPS_INT>(m_ordering);
        instance.icntl[15] = static_cast<MUMPS_INT>(m_maxThreads);
        instance.icntl[13] = static_cast<MUMPS_INT>(m_workspacePercentage);
      }

      /**
       * @brief Records MUMPS's status and the stage that survived it.
       *
       * A failing analysis, factorization or solve is an outcome rather than a
       * defect in the caller, so it is reported through getInfo() instead of
       * raising.
       */
      void record(Optional<Factorization> factorization)
      {
        const auto& instance = m_resources.instance;
        m_info.status = static_cast<Integer>(instance.infog[0]);
        m_info.success = instance.infog[0] >= 0;
        m_info.factorization = factorization;
      }

      void check(StringView operation) const
      {
        const auto& instance = m_resources.instance;
        if (instance.infog[0] < 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MUMPS " << operation
            << " failed with INFOG(1) = " << static_cast<Integer>(instance.infog[0])
            << " and INFOG(2) = " << static_cast<Integer>(instance.infog[1]) << "."
            << Alert::Raise;
        }
      }

      /// @brief Outcome of the most recent operation.
      Info m_info;

      Index m_maxThreads = 0;
      Symmetry m_symmetry = Symmetry::Unsymmetric;
      Ordering m_ordering = Ordering::Automatic;
      Integer m_workspacePercentage = 0;
      Array<MUMPS_INT> m_rowIndices;
      Array<MUMPS_INT> m_columnIndices;
      Math::Vector<ScalarType> m_values;
      Resources m_resources;
  };
}

#endif
#endif
