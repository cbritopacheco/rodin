/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_BOOMERAMG_H
#define RODIN_KOKKOS_SOLVER_BOOMERAMG_H

/**
 * @file
 * @brief Conjugate gradient of hypre preconditioned by BoomerAMG, for a
 * Kokkos::Math::SparseMatrix.
 */

#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include <HYPRE.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Solver/LinearSolver.h"

#include "Rodin/Kokkos/Math/LinearSystem.h"

namespace Rodin::Solver
{
  /**
   * @brief Conjugate gradient of hypre preconditioned by one BoomerAMG
   * V-cycle.
   * @tparam LinearSystem Type of linear system to solve.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref BoomerAMG "BoomerAMG<::Rodin::Math::LinearSystem<Kokkos::Math::SparseMatrix<Real, Backend>, ::Rodin::Math::Vector<Real>>>" | hypre on the arrays of a device matrix. |
   *
   * @see @ref BoomerAMGSpecializations "BoomerAMG Template Specializations"
   */
  template <class LinearSystem>
  class BoomerAMG;

  /**
   * @defgroup BoomerAMGSpecializations BoomerAMG Template Specializations
   * @brief Template specializations of the BoomerAMG class.
   * @see @ref BoomerAMG
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for BoomerAMG.
   * @param problem Problem whose linear system the solver solves.
   */
  template <class LinearSystem>
  BoomerAMG(Variational::ProblemBase<LinearSystem>& problem) -> BoomerAMG<LinearSystem>;

  /**
   * @ingroup BoomerAMGSpecializations
   * @brief Conjugate gradient of hypre preconditioned by BoomerAMG, for a
   * system whose matrix is a Kokkos::Math::SparseMatrix.
   *
   * The matrix must be symmetric positive definite. hypre copies the CSR
   * arrays of the device matrix in the memory of the backend, sets up the
   * hierarchy and solves there; only the right-hand side and the solution
   * cross to the host. A solve performs at least one iteration, even from the
   * solution.
   *
   * The hierarchy is retained from one solve to the next, together with the
   * copy of the matrix it was built from, as PETSc does with
   * @c KSPSetReusePreconditioner. The ordinary @ref solve entry point builds
   * a hierarchy only when none is retained; the conjugate gradient always
   * runs on the current matrix, so a hierarchy built from an older matrix of
   * the same size still gives the solution of the current system, in more
   * iterations as the matrix drifts away. Call @ref setup to build the
   * hierarchy again from the current matrix, and @ref clear to release it.
   * Changing the coarsening settings releases the hierarchy.
   *
   * The solve stops when @f$ \|b - Ax\| \leq \text{tol} \|b\| @f$, with a
   * tolerance of @f$ 10^{-10} @f$ and at most 1000 iterations by default.
   * The hierarchy uses settings that run on a GPU: PMIS coarsening, extended
   * plus i interpolation, one sweep of l1 Jacobi and Gaussian elimination on
   * the coarsest grid. The strength threshold (0.25 by default), the levels
   * of aggressive coarsening (none by default) and the interpolation
   * truncation (4 elements per row by default) are adjustable; three
   * dimensional problems usually need a threshold of 0.5 to 0.7 and
   * aggressive coarsening to fit in the memory of the device.
   *
   * MPI and hypre are initialized on the first setup when the application
   * has not done so. Their finalization is left to the application, after the
   * last solver is destroyed, since a solver releases its hierarchy on
   * destruction. The solve runs on @c MPI_COMM_SELF.
   *
   * @tparam Backend Kokkos backend of the matrix.
   */
  template <class Backend, class VectorType_>
  class BoomerAMG<::Rodin::Math::LinearSystem<
    ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, VectorType_>>
    final : public LinearSolverBase<::Rodin::Math::LinearSystem<
              ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, VectorType_>>
  {
    static_assert(std::is_same_v<HYPRE_Int, ::Rodin::Math::SparseIndex> &&
        std::is_same_v<HYPRE_BigInt, ::Rodin::Math::SparseIndex>,
      "hypre must be built with indices of the size of ::Rodin::Math::SparseIndex.");
    static_assert(
      std::is_same_v<HYPRE_Complex, Real>, "hypre must be built in double precision.");

    public:
      /// @brief Scalar type.
      using ScalarType = Real;

      /// @brief Host vector type.
      using VectorType = VectorType_;

      /// @brief Device matrix type.
      using OperatorType = ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>;

      /// @brief Linear system type.
      using LinearSystemType = ::Rodin::Math::LinearSystem<OperatorType, VectorType>;

      /// @brief Problem type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystemType>;

      /// @brief Execution space of the matrix.
      using ExecutionSpace = typename OperatorType::ExecutionSpace;

      /// @brief Memory space of the matrix.
      using MemorySpace = typename OperatorType::MemorySpace;

      using Parent::solve;

      /**
       * @brief Constructs the solver of a problem.
       * @param problem Problem to solve.
       */
      explicit BoomerAMG(ProblemBaseType& problem)
        : Parent(problem)
      {}

      /**
       * @brief Copies the settings and the outcome of the last solve, without
       * the hierarchy.
       * @param other Solver to copy.
       */
      BoomerAMG(const BoomerAMG& other)
        : Parent(other),
          m_tolerance(other.m_tolerance),
          m_maxIterations(other.m_maxIterations),
          m_strongThreshold(other.m_strongThreshold),
          m_aggressiveLevels(other.m_aggressiveLevels),
          m_maxInterpolationElements(other.m_maxInterpolationElements),
          m_converged(other.m_converged),
          m_iterations(other.m_iterations),
          m_error(other.m_error)
      {}

      /**
       * @brief Moves the solver.
       * @param other Solver to move.
       */
      BoomerAMG(BoomerAMG&& other) = default;

      /// @brief Destructor.
      ~BoomerAMG() override = default;

      /**
       * @brief Sets the relative tolerance.
       * @param tolerance Tolerance on @f$ \|b - Ax\| / \|b\| @f$.
       * @returns This solver.
       */
      BoomerAMG& setTolerance(Real tolerance)
      {
        m_tolerance = tolerance;
        return *this;
      }

      /**
       * @brief Sets the iteration limit of the conjugate gradient.
       * @param count Maximum number of iterations.
       * @returns This solver.
       */
      BoomerAMG& setMaxIterations(size_t count)
      {
        m_maxIterations = count;
        return *this;
      }

      /**
       * @brief Sets the strength threshold of the coarsening.
       * @param threshold Threshold in @f$ (0, 1) @f$.
       * Releases the retained hierarchy.
       * @returns This solver.
       */
      BoomerAMG& setStrongThreshold(Real threshold)
      {
        m_strongThreshold = threshold;
        clear();
        return *this;
      }

      /**
       * @brief Sets the number of levels of aggressive coarsening.
       * @param count Number of levels, from the finest.
       * Releases the retained hierarchy.
       * @returns This solver.
       */
      BoomerAMG& setAggressiveLevels(size_t count)
      {
        m_aggressiveLevels = count;
        clear();
        return *this;
      }

      /**
       * @brief Sets the maximum number of elements per row of the
       * interpolation.
       * @param count Number of elements, or zero for no truncation.
       * Releases the retained hierarchy.
       * @returns This solver.
       */
      BoomerAMG& setMaxInterpolationElements(size_t count)
      {
        m_maxInterpolationElements = count;
        clear();
        return *this;
      }

      /**
       * @brief Builds the hierarchy from the current matrix of a system.
       * @param system System whose matrix is copied for the hierarchy.
       *
       * Replaces any retained hierarchy.
       */
      void setup(const LinearSystemType& system)
      {
        const auto& matrix = system.getOperator();
        validate(matrix);
        clear();
        const size_t n = matrix.getRowCount();
        if (n == 0)
          return;
        initializeSession();
        configure();
        // The setup only reads the partition of the vectors.
        const Array<Real> zero("Rodin::Solver::BoomerAMG::zero", n);
        ::Kokkos::fence();
        Objects objects;
        const HYPRE_ParVector B = createVector(objects.b, zero);
        const HYPRE_ParVector X = createVector(objects.x, zero);
        build(matrix, B, X);
      }

      /// @brief Releases the retained hierarchy.
      void clear() noexcept
      {
        m_hierarchy.reset();
      }

      /**
       * @brief Whether a hierarchy is retained.
       * @returns True if a hierarchy has been set up and not cleared.
       */
      Boolean isSetup() const noexcept
      {
        return static_cast<bool>(m_hierarchy);
      }

      /**
       * @brief Solves the system.
       * @param[in,out] system System, whose solution is the initial guess when
       * its size matches.
       *
       * Builds the hierarchy from the current matrix when none is retained.
       */
      void solve(LinearSystemType& system) override
      {
        const auto& matrix = system.getOperator();
        validate(matrix);
        const size_t n = matrix.getRowCount();
        if (static_cast<size_t>(system.getVector().size()) != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }
        if (m_hierarchy && m_hierarchy->rows != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The matrix size changed since the setup of the hierarchy."
            << Alert::Raise;
        }
        // hypre returns a zero solution for a zero right-hand side without
        // reporting convergence.
        if (n == 0 || system.getVector().isZero(0))
        {
          m_converged = true;
          m_iterations = 0;
          m_error = 0;
          system.getSolution().setZero(static_cast<Eigen::Index>(n));
          return;
        }
        initializeSession();
        configure();

        Array<Real> b(::Kokkos::view_alloc(
                        ::Kokkos::WithoutInitializing, "Rodin::Solver::BoomerAMG::b"),
          n);
        Array<Real> x(::Kokkos::view_alloc(
                        ::Kokkos::WithoutInitializing, "Rodin::Solver::BoomerAMG::x"),
          n);
        if constexpr (requires { system.getSolution().getView(); })
        {
          if (static_cast<size_t>(system.getSolution().size()) != n)
            system.getSolution().setZero(n);
        }
        upload(system.getVector(), b);
        if (static_cast<size_t>(system.getSolution().size()) == n)
          upload(system.getSolution(), x);
        else
          ::Kokkos::deep_copy(x, Real(0));
        // hypre runs on its own streams.
        ::Kokkos::fence();

        Objects objects;
        const HYPRE_ParVector B = createVector(objects.b, b);
        const HYPRE_ParVector X = createVector(objects.x, x);
        // A hierarchy built here comes from the current matrix, whose copy
        // serves the conjugate gradient as well.
        HYPRE_ParCSRMatrix A = nullptr;
        if (m_hierarchy)
          A = createMatrix(objects.matrix, matrix);
        else
          A = build(matrix, B, X);

        check(
          HYPRE_ParCSRPCGCreate(MPI_COMM_SELF, &objects.pcg), "HYPRE_ParCSRPCGCreate");
        check(HYPRE_PCGSetTol(objects.pcg, m_tolerance), "HYPRE_PCGSetTol");
        check(HYPRE_PCGSetMaxIter(objects.pcg, static_cast<HYPRE_Int>(m_maxIterations)),
          "HYPRE_PCGSetMaxIter");
        check(HYPRE_PCGSetTwoNorm(objects.pcg, 1), "HYPRE_PCGSetTwoNorm");
        // The V-cycle runs on the matrix of the hierarchy, which the setup of
        // the conjugate gradient leaves untouched.
        check(HYPRE_PCGSetPrecond(objects.pcg,
                reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSolve), &skipSetup,
                m_hierarchy->amg),
          "HYPRE_PCGSetPrecond");
        check(HYPRE_PCGSetPrecondMatrix(
                objects.pcg, reinterpret_cast<HYPRE_Matrix>(m_hierarchy->A)),
          "HYPRE_PCGSetPrecondMatrix");
        check(HYPRE_ParCSRPCGSetup(objects.pcg, A, B, X), "HYPRE_ParCSRPCGSetup");
        // A solve that does not converge reports it through the error flag,
        // which is an outcome rather than a defect.
        const HYPRE_Int status = HYPRE_ParCSRPCGSolve(objects.pcg, A, B, X);
        check(status & ~HYPRE_ERROR_CONV, "HYPRE_ParCSRPCGSolve");
        HYPRE_ClearAllErrors();

        HYPRE_Int iterations = 0;
        HYPRE_Int converged = 0;
        HYPRE_Real error = 0;
        check(HYPRE_PCGGetNumIterations(objects.pcg, &iterations),
          "HYPRE_PCGGetNumIterations");
        check(HYPRE_PCGGetConverged(objects.pcg, &converged), "HYPRE_PCGGetConverged");
        check(HYPRE_PCGGetFinalRelativeResidualNorm(objects.pcg, &error),
          "HYPRE_PCGGetFinalRelativeResidualNorm");
        m_iterations = static_cast<size_t>(iterations);
        m_converged = converged != 0;
        m_error = static_cast<Real>(error);

        check(HYPRE_IJVectorGetValues(
                objects.x, static_cast<HYPRE_Int>(n), nullptr, x.data()),
          "HYPRE_IJVectorGetValues");
        ::Kokkos::fence();
        download(x, system.getSolution());
      }

      /**
       * @brief Whether the last solve reached the tolerance.
       * @returns True if the last solve converged.
       */
      Boolean success() const
      {
        return m_converged;
      }

      /**
       * @brief Number of iterations of the last solve.
       * @returns Iteration count reported by hypre.
       */
      size_t getIterationNumber() const
      {
        return m_iterations;
      }

      /**
       * @brief Relative residual of the last solve.
       * @returns Final relative residual norm reported by hypre.
       */
      Real getError() const
      {
        return m_error;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new BoomerAMG instance.
       */
      BoomerAMG* copy() const noexcept override
      {
        return new BoomerAMG(*this);
      }

      /// @brief Writes the row count and the row number of each row of a CSR matrix.
      template <class Starts, class Counts, class Rows>
      struct FillRows
      {
        /// @brief CSR row starts.
        Starts starts;

        /// @brief Row counts, one per row.
        Counts counts;

        /// @brief Global row numbers, one per row.
        Rows rows;

        /**
         * @brief Writes the count and the number of one row.
         * @param i Row index.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          counts(i) = starts(i + 1) - starts(i);
          rows(i) = static_cast<HYPRE_BigInt>(i);
        }
      };

      /**
       * @brief Writes the row counts and the row numbers of a CSR matrix.
       * @param starts Row starts.
       * @param counts Row counts, one per row.
       * @param rows Row numbers, one per row.
       */
      template <class Starts, class Counts, class Rows>
      static void fillRows(const Starts& starts, const Counts& counts, const Rows& rows)
      {
        ::Kokkos::parallel_for("Rodin::Solver::BoomerAMG::rows",
          ::Kokkos::RangePolicy<ExecutionSpace>(0, counts.extent(0)),
          FillRows<Starts, Counts, Rows>{starts, counts, rows});
      }

    private:
      template <class T>
      using Array = ::Kokkos::View<T*, MemorySpace>;

      /// @brief Memory of hypre objects that read the arrays of the backend.
      static constexpr HYPRE_MemoryLocation s_memoryLocation =
        ::Kokkos::SpaceAccessibility<::Kokkos::HostSpace, MemorySpace>::accessible
        ? HYPRE_MEMORY_HOST
        : HYPRE_MEMORY_DEVICE;

      /// @brief hypre objects of one solve, destroyed in reverse order.
      struct Objects
      {
        HYPRE_IJMatrix matrix = nullptr;
        HYPRE_IJVector b = nullptr;
        HYPRE_IJVector x = nullptr;
        HYPRE_Solver pcg = nullptr;

        Objects() = default;

        /**
         * @brief Not copyable, the objects being destroyed once.
         * @param other Objects not copied.
         */
        Objects(const Objects& other) = delete;

        /**
         * @brief Not copyable, the objects being destroyed once.
         * @param other Objects not copied.
         * @returns Reference to these objects.
         */
        Objects& operator=(const Objects& other) = delete;

        ~Objects()
        {
          if (pcg)
            HYPRE_ParCSRPCGDestroy(pcg);
          if (x)
            HYPRE_IJVectorDestroy(x);
          if (b)
            HYPRE_IJVectorDestroy(b);
          if (matrix)
            HYPRE_IJMatrixDestroy(matrix);
        }
      };

      /// @brief Retained hierarchy and the copy of the matrix it comes from.
      struct Hierarchy
      {
        HYPRE_IJMatrix matrix = nullptr;
        HYPRE_ParCSRMatrix A = nullptr;
        HYPRE_Solver amg = nullptr;
        size_t rows = 0;

        Hierarchy() = default;

        /**
         * @brief Not copyable, the hierarchy being destroyed once.
         * @param other Hierarchy not copied.
         */
        Hierarchy(const Hierarchy& other) = delete;

        /**
         * @brief Not copyable, the hierarchy being destroyed once.
         * @param other Hierarchy not copied.
         * @returns Reference to this hierarchy.
         */
        Hierarchy& operator=(const Hierarchy& other) = delete;

        ~Hierarchy()
        {
          if (amg)
            HYPRE_BoomerAMGDestroy(amg);
          if (matrix)
            HYPRE_IJMatrixDestroy(matrix);
        }
      };

      /**
       * @brief Preconditioner setup that the conjugate gradient calls, the
       * hierarchy being built beforehand.
       * @param solver Preconditioner, already set up.
       * @param A Matrix of the system.
       * @param b Right-hand side.
       * @param x Solution.
       * @returns Zero, the success code of hypre.
       */
      static HYPRE_Int skipSetup([[maybe_unused]] HYPRE_Solver solver,
        [[maybe_unused]] HYPRE_Matrix A, [[maybe_unused]] HYPRE_Vector b,
        [[maybe_unused]] HYPRE_Vector x)
      {
        return 0;
      }

      /**
       * @brief Checks that hypre can take a matrix.
       * @param matrix Matrix of the system.
       */
      void validate(const OperatorType& matrix) const
      {
        const size_t n = matrix.getRowCount();
        if (matrix.getColumnCount() != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "BoomerAMG requires a square matrix." << Alert::Raise;
        }
        if (n > static_cast<size_t>(std::numeric_limits<HYPRE_Int>::max()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The system is too large for the indices of hypre." << Alert::Raise;
        }
      }

      /// @brief Places the next hypre objects and computations on the backend.
      void configure() const
      {
        check(HYPRE_SetMemoryLocation(s_memoryLocation), "HYPRE_SetMemoryLocation");
        check(
          HYPRE_SetExecutionPolicy(
            s_memoryLocation == HYPRE_MEMORY_HOST ? HYPRE_EXEC_HOST : HYPRE_EXEC_DEVICE),
          "HYPRE_SetExecutionPolicy");
      }

      /**
       * @brief Creates a hypre matrix holding a copy of the device matrix.
       * @param ij Created matrix.
       * @param matrix Device matrix.
       * @returns Parallel matrix of @p ij.
       */
      HYPRE_ParCSRMatrix createMatrix(
        HYPRE_IJMatrix& ij, const OperatorType& matrix) const
      {
        // Row counts and row numbers, which the IJ interface takes with the
        // row starts.
        const size_t n = matrix.getRowCount();
        const auto rowCount = static_cast<HYPRE_Int>(n);
        Array<HYPRE_Int> counts(::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                                  "Rodin::Solver::BoomerAMG::counts"),
          n);
        Array<HYPRE_BigInt> rows(::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                                   "Rodin::Solver::BoomerAMG::rows"),
          n);
        const auto starts = matrix.getRowStarts();
        fillRows(starts, counts, rows);
        ::Kokkos::fence();

        const HYPRE_BigInt last = rowCount - 1;
        check(HYPRE_IJMatrixCreate(MPI_COMM_SELF, 0, last, 0, last, &ij),
          "HYPRE_IJMatrixCreate");
        check(
          HYPRE_IJMatrixSetObjectType(ij, HYPRE_PARCSR), "HYPRE_IJMatrixSetObjectType");
        check(HYPRE_IJMatrixInitialize_v2(ij, s_memoryLocation),
          "HYPRE_IJMatrixInitialize_v2");
        check(HYPRE_IJMatrixSetValues2(ij, rowCount, counts.data(), rows.data(),
                starts.data(), matrix.getColumns().data(), matrix.getValues().data()),
          "HYPRE_IJMatrixSetValues2");
        check(HYPRE_IJMatrixAssemble(ij), "HYPRE_IJMatrixAssemble");
        // The row counts are released on return.
        ::Kokkos::fence();
        HYPRE_ParCSRMatrix res = nullptr;
        check(HYPRE_IJMatrixGetObject(ij, reinterpret_cast<void**>(&res)),
          "HYPRE_IJMatrixGetObject");
        return res;
      }

      /**
       * @brief Builds and retains the hierarchy of a matrix.
       * @param matrix Device matrix.
       * @param B Right-hand side, whose partition the setup reads.
       * @param X Solution, whose partition the setup reads.
       * @returns Copy of @p matrix in the hierarchy.
       */
      HYPRE_ParCSRMatrix build(
        const OperatorType& matrix, HYPRE_ParVector B, HYPRE_ParVector X)
      {
        auto hierarchy = std::make_unique<Hierarchy>();
        hierarchy->rows = matrix.getRowCount();
        hierarchy->A = createMatrix(hierarchy->matrix, matrix);
        HYPRE_Solver& amg = hierarchy->amg;
        check(HYPRE_BoomerAMGCreate(&amg), "HYPRE_BoomerAMGCreate");
        check(HYPRE_BoomerAMGSetCoarsenType(amg, 8), "HYPRE_BoomerAMGSetCoarsenType");
        check(HYPRE_BoomerAMGSetInterpType(amg, 6), "HYPRE_BoomerAMGSetInterpType");
        check(HYPRE_BoomerAMGSetRelaxType(amg, 18), "HYPRE_BoomerAMGSetRelaxType");
        check(HYPRE_BoomerAMGSetCycleRelaxType(amg, 9, 3),
          "HYPRE_BoomerAMGSetCycleRelaxType");
        check(HYPRE_BoomerAMGSetNumSweeps(amg, 1), "HYPRE_BoomerAMGSetNumSweeps");
        check(HYPRE_BoomerAMGSetStrongThreshold(amg, m_strongThreshold),
          "HYPRE_BoomerAMGSetStrongThreshold");
        check(
          HYPRE_BoomerAMGSetAggNumLevels(amg, static_cast<HYPRE_Int>(m_aggressiveLevels)),
          "HYPRE_BoomerAMGSetAggNumLevels");
        check(HYPRE_BoomerAMGSetPMaxElmts(
                amg, static_cast<HYPRE_Int>(m_maxInterpolationElements)),
          "HYPRE_BoomerAMGSetPMaxElmts");
        check(HYPRE_BoomerAMGSetMaxIter(amg, 1), "HYPRE_BoomerAMGSetMaxIter");
        check(HYPRE_BoomerAMGSetTol(amg, 0), "HYPRE_BoomerAMGSetTol");
        check(HYPRE_BoomerAMGSetup(amg, hierarchy->A, B, X), "HYPRE_BoomerAMGSetup");
        m_hierarchy = std::move(hierarchy);
        return m_hierarchy->A;
      }

      /**
       * @brief Raises the error of a hypre call.
       * @param error Error flag returned by the call.
       * @param call Name of the call.
       */
      void check(HYPRE_Int error, const char* call) const
      {
        if (error == 0)
          return;
        char description[256] = {};
        HYPRE_DescribeError(error, description);
        HYPRE_ClearAllErrors();
        Alert::MemberFunctionException(*this, __func__)
          << call << " failed: " << description << Alert::Raise;
      }

      /// @brief Ensures that MPI and hypre are initialized.
      void initializeSession() const
      {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (finalized || HYPRE_Finalized())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MPI or hypre has been finalized, so hypre can no longer be used."
            << Alert::Raise;
        }
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized)
        {
          int argc = 0;
          char** argv = nullptr;
          int provided = 0;
          if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) !=
            MPI_SUCCESS)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "MPI initialization failed." << Alert::Raise;
          }
        }
        if (!HYPRE_Initialized())
          check(HYPRE_Initialize(), "HYPRE_Initialize");
      }

      /**
       * @brief Creates a hypre vector holding a copy of an array.
       * @param vector Created vector.
       * @param values Values in the memory of the backend.
       * @returns Parallel vector of @p vector.
       */
      HYPRE_ParVector createVector(
        HYPRE_IJVector& vector, const Array<Real>& values) const
      {
        const auto last = static_cast<HYPRE_BigInt>(values.extent(0)) - 1;
        check(
          HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, last, &vector), "HYPRE_IJVectorCreate");
        check(HYPRE_IJVectorSetObjectType(vector, HYPRE_PARCSR),
          "HYPRE_IJVectorSetObjectType");
        check(HYPRE_IJVectorInitialize_v2(vector, s_memoryLocation),
          "HYPRE_IJVectorInitialize_v2");
        check(HYPRE_IJVectorSetValues(
                vector, static_cast<HYPRE_Int>(values.extent(0)), nullptr, values.data()),
          "HYPRE_IJVectorSetValues");
        check(HYPRE_IJVectorAssemble(vector), "HYPRE_IJVectorAssemble");
        HYPRE_ParVector res = nullptr;
        check(HYPRE_IJVectorGetObject(vector, reinterpret_cast<void**>(&res)),
          "HYPRE_IJVectorGetObject");
        return res;
      }

      /// @param source Backend coefficients. @param destination Borrowed backend view.
      static void upload(const ::Rodin::Kokkos::Math::Vector<Real, Backend>& source,
        Array<Real>& destination)
      {
        destination = source.getView();
      }

      /// @param source Backend iterate. @param destination Backend solution owner.
      static void download(const Array<Real>& source,
        ::Rodin::Kokkos::Math::Vector<Real, Backend>& destination)
      {
        if (source.data() == destination.data())
          return;
        destination.resize(source.extent(0));
        ::Kokkos::deep_copy(destination.getView(), source);
      }

      /**
       * @brief Copies a host vector to the memory of the backend.
       * @param host Host vector.
       * @param device Array of the same size.
       */
      static void upload(
        const ::Rodin::Math::Vector<Real>& host, const Array<Real>& device)
      {
        const ::Kokkos::View<const Real*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged>
          h(host.data(), static_cast<size_t>(host.size()));
        ::Kokkos::deep_copy(device, h);
      }

      /**
       * @brief Copies an array of the backend to the host.
       * @param device Array.
       * @param host Host vector, resized.
       */
      static void download(const Array<Real>& device, ::Rodin::Math::Vector<Real>& host)
      {
        host.resize(static_cast<Eigen::Index>(device.extent(0)));
        const ::Kokkos::View<Real*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> h(
          host.data(), device.extent(0));
        ::Kokkos::deep_copy(h, device);
      }

      Real m_tolerance = 1e-10;
      size_t m_maxIterations = 1000;
      Real m_strongThreshold = 0.25;
      size_t m_aggressiveLevels = 0;
      size_t m_maxInterpolationElements = 4;
      Boolean m_converged = false;
      size_t m_iterations = 0;
      Real m_error = 0;
      std::unique_ptr<Hierarchy> m_hierarchy;
  };
}

#endif
