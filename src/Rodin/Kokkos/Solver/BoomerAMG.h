/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_BOOMERAMG_H
#define RODIN_KOKKOS_SOLVER_BOOMERAMG_H

/**
 * @file
 * @brief Conjugate gradient of hypre preconditioned by BoomerAMG, for a
 * Kokkos::SparseMatrix.
 */

#include <limits>
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

#include "Rodin/Kokkos/LinearSystem.h"

namespace Rodin::Solver
{
  /**
   * @brief Conjugate gradient of hypre preconditioned by one BoomerAMG
   * V-cycle.
   * @tparam LinearSystem Type of linear system to solve.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref BoomerAMG "BoomerAMG<Math::LinearSystem<Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>>" | hypre on the arrays of a device matrix. |
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
   */
  template <class LinearSystem>
  BoomerAMG(Variational::ProblemBase<LinearSystem>&) -> BoomerAMG<LinearSystem>;

  /**
   * @ingroup BoomerAMGSpecializations
   * @brief Conjugate gradient of hypre preconditioned by BoomerAMG, for a
   * system whose matrix is a Kokkos::SparseMatrix.
   *
   * The matrix must be symmetric positive definite. hypre copies the CSR
   * arrays of the device matrix in the memory of the backend, sets up the
   * hierarchy and solves there; only the right-hand side and the solution
   * cross to the host. Every solve sets the hierarchy up again and performs at
   * least one iteration, even from the solution.
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
   * MPI and hypre are initialized on the first solve when the application
   * has not done so. Their finalization is left to the application, after the
   * last solver is destroyed. The solve runs on @c MPI_COMM_SELF.
   *
   * @tparam Backend Kokkos backend of the matrix.
   */
  template <class Backend>
  class BoomerAMG<
    Math::LinearSystem<::Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>>
    final : public LinearSolverBase<Math::LinearSystem<
              ::Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>>
  {
    static_assert(std::is_same_v<HYPRE_Int, Math::SparseIndex> &&
        std::is_same_v<HYPRE_BigInt, Math::SparseIndex>,
      "hypre must be built with indices of the size of Math::SparseIndex.");
    static_assert(
      std::is_same_v<HYPRE_Complex, Real>, "hypre must be built in double precision.");

    public:
      /// @brief Scalar type.
      using ScalarType = Real;

      /// @brief Host vector type.
      using VectorType = Math::Vector<Real>;

      /// @brief Device matrix type.
      using OperatorType = ::Rodin::Kokkos::SparseMatrix<Real, Backend>;

      /// @brief Linear system type.
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

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
       * @brief Copies the settings and the outcome of the last solve.
       * @param other Solver to copy.
       */
      BoomerAMG(const BoomerAMG& other) = default;

      /// @brief Moves the solver.
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
       * @returns This solver.
       */
      BoomerAMG& setStrongThreshold(Real threshold)
      {
        m_strongThreshold = threshold;
        return *this;
      }

      /**
       * @brief Sets the number of levels of aggressive coarsening.
       * @param count Number of levels, from the finest.
       * @returns This solver.
       */
      BoomerAMG& setAggressiveLevels(size_t count)
      {
        m_aggressiveLevels = count;
        return *this;
      }

      /**
       * @brief Sets the maximum number of elements per row of the
       * interpolation.
       * @param count Number of elements, or zero for no truncation.
       * @returns This solver.
       */
      BoomerAMG& setMaxInterpolationElements(size_t count)
      {
        m_maxInterpolationElements = count;
        return *this;
      }

      /**
       * @brief Solves the system.
       * @param[in,out] system System, whose solution is the initial guess when
       * its size matches.
       */
      void solve(LinearSystemType& system) override
      {
        const auto& matrix = system.getOperator();
        const size_t n = matrix.getRowCount();
        if (matrix.getColumnCount() != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "BoomerAMG requires a square matrix." << Alert::Raise;
        }
        if (static_cast<size_t>(system.getVector().size()) != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }
        if (n > static_cast<size_t>(std::numeric_limits<HYPRE_Int>::max()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The system is too large for the indices of hypre." << Alert::Raise;
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

        const auto location = s_memoryLocation;
        check(HYPRE_SetMemoryLocation(location), "HYPRE_SetMemoryLocation");
        check(HYPRE_SetExecutionPolicy(
                location == HYPRE_MEMORY_HOST ? HYPRE_EXEC_HOST : HYPRE_EXEC_DEVICE),
          "HYPRE_SetExecutionPolicy");

        // Row counts and row numbers, which the IJ interface takes with the
        // row starts.
        const auto rowCount = static_cast<HYPRE_Int>(n);
        Array<HYPRE_Int> counts(::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                                  "Rodin::Solver::BoomerAMG::counts"),
          n);
        Array<HYPRE_BigInt> rows(::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                                   "Rodin::Solver::BoomerAMG::rows"),
          n);
        const auto starts = matrix.getRowStarts();
        ::Kokkos::parallel_for(
          "Rodin::Solver::BoomerAMG::rows", ::Kokkos::RangePolicy<ExecutionSpace>(0, n),
          KOKKOS_LAMBDA(const size_t i) {
            counts(i) = starts(i + 1) - starts(i);
            rows(i) = static_cast<HYPRE_BigInt>(i);
          });

        Array<Real> b(::Kokkos::view_alloc(
                        ::Kokkos::WithoutInitializing, "Rodin::Solver::BoomerAMG::b"),
          n);
        Array<Real> x(::Kokkos::view_alloc(
                        ::Kokkos::WithoutInitializing, "Rodin::Solver::BoomerAMG::x"),
          n);
        upload(system.getVector(), b);
        if (static_cast<size_t>(system.getSolution().size()) == n)
          upload(system.getSolution(), x);
        else
          ::Kokkos::deep_copy(x, Real(0));
        // hypre runs on its own streams.
        ::Kokkos::fence();

        Objects objects;
        const HYPRE_BigInt last = rowCount - 1;
        check(HYPRE_IJMatrixCreate(MPI_COMM_SELF, 0, last, 0, last, &objects.matrix),
          "HYPRE_IJMatrixCreate");
        check(HYPRE_IJMatrixSetObjectType(objects.matrix, HYPRE_PARCSR),
          "HYPRE_IJMatrixSetObjectType");
        check(HYPRE_IJMatrixInitialize_v2(objects.matrix, location),
          "HYPRE_IJMatrixInitialize_v2");
        check(
          HYPRE_IJMatrixSetValues2(objects.matrix, rowCount, counts.data(), rows.data(),
            starts.data(), matrix.getColumns().data(), matrix.getValues().data()),
          "HYPRE_IJMatrixSetValues2");
        check(HYPRE_IJMatrixAssemble(objects.matrix), "HYPRE_IJMatrixAssemble");
        HYPRE_ParCSRMatrix A = nullptr;
        check(HYPRE_IJMatrixGetObject(objects.matrix, reinterpret_cast<void**>(&A)),
          "HYPRE_IJMatrixGetObject");

        const HYPRE_ParVector B = createVector(objects.b, b, last, location);
        const HYPRE_ParVector X = createVector(objects.x, x, last, location);

        check(HYPRE_BoomerAMGCreate(&objects.amg), "HYPRE_BoomerAMGCreate");
        check(
          HYPRE_BoomerAMGSetCoarsenType(objects.amg, 8), "HYPRE_BoomerAMGSetCoarsenType");
        check(
          HYPRE_BoomerAMGSetInterpType(objects.amg, 6), "HYPRE_BoomerAMGSetInterpType");
        check(
          HYPRE_BoomerAMGSetRelaxType(objects.amg, 18), "HYPRE_BoomerAMGSetRelaxType");
        check(HYPRE_BoomerAMGSetCycleRelaxType(objects.amg, 9, 3),
          "HYPRE_BoomerAMGSetCycleRelaxType");
        check(HYPRE_BoomerAMGSetNumSweeps(objects.amg, 1), "HYPRE_BoomerAMGSetNumSweeps");
        check(HYPRE_BoomerAMGSetStrongThreshold(objects.amg, m_strongThreshold),
          "HYPRE_BoomerAMGSetStrongThreshold");
        check(HYPRE_BoomerAMGSetAggNumLevels(
                objects.amg, static_cast<HYPRE_Int>(m_aggressiveLevels)),
          "HYPRE_BoomerAMGSetAggNumLevels");
        check(HYPRE_BoomerAMGSetPMaxElmts(
                objects.amg, static_cast<HYPRE_Int>(m_maxInterpolationElements)),
          "HYPRE_BoomerAMGSetPMaxElmts");
        check(HYPRE_BoomerAMGSetMaxIter(objects.amg, 1), "HYPRE_BoomerAMGSetMaxIter");
        check(HYPRE_BoomerAMGSetTol(objects.amg, 0), "HYPRE_BoomerAMGSetTol");

        check(
          HYPRE_ParCSRPCGCreate(MPI_COMM_SELF, &objects.pcg), "HYPRE_ParCSRPCGCreate");
        check(HYPRE_PCGSetTol(objects.pcg, m_tolerance), "HYPRE_PCGSetTol");
        check(HYPRE_PCGSetMaxIter(objects.pcg, static_cast<HYPRE_Int>(m_maxIterations)),
          "HYPRE_PCGSetMaxIter");
        check(HYPRE_PCGSetTwoNorm(objects.pcg, 1), "HYPRE_PCGSetTwoNorm");
        check(
          HYPRE_PCGSetPrecond(objects.pcg,
            reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSolve),
            reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSetup), objects.amg),
          "HYPRE_PCGSetPrecond");
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

        check(HYPRE_IJVectorGetValues(objects.x, rowCount, nullptr, x.data()),
          "HYPRE_IJVectorGetValues");
        ::Kokkos::fence();
        download(x, system.getSolution());
      }

      /// @brief Whether the last solve reached the tolerance.
      Boolean success() const
      {
        return m_converged;
      }

      /// @brief Number of iterations of the last solve.
      size_t getIterationNumber() const
      {
        return m_iterations;
      }

      /// @brief Relative residual of the last solve.
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
        HYPRE_Solver amg = nullptr;
        HYPRE_Solver pcg = nullptr;

        Objects() = default;

        Objects(const Objects&) = delete;

        Objects& operator=(const Objects&) = delete;

        ~Objects()
        {
          if (pcg)
            HYPRE_ParCSRPCGDestroy(pcg);
          if (amg)
            HYPRE_BoomerAMGDestroy(amg);
          if (x)
            HYPRE_IJVectorDestroy(x);
          if (b)
            HYPRE_IJVectorDestroy(b);
          if (matrix)
            HYPRE_IJMatrixDestroy(matrix);
        }
      };

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

      /**
       * @brief Ensures that MPI and hypre are initialized.
       */
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
       * @param last Last row.
       * @param location Memory of the vector.
       * @returns Parallel vector of @p vector.
       */
      HYPRE_ParVector createVector(HYPRE_IJVector& vector, const Array<Real>& values,
        HYPRE_BigInt last, HYPRE_MemoryLocation location) const
      {
        check(
          HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, last, &vector), "HYPRE_IJVectorCreate");
        check(HYPRE_IJVectorSetObjectType(vector, HYPRE_PARCSR),
          "HYPRE_IJVectorSetObjectType");
        check(
          HYPRE_IJVectorInitialize_v2(vector, location), "HYPRE_IJVectorInitialize_v2");
        check(HYPRE_IJVectorSetValues(
                vector, static_cast<HYPRE_Int>(values.extent(0)), nullptr, values.data()),
          "HYPRE_IJVectorSetValues");
        check(HYPRE_IJVectorAssemble(vector), "HYPRE_IJVectorAssemble");
        HYPRE_ParVector res = nullptr;
        check(HYPRE_IJVectorGetObject(vector, reinterpret_cast<void**>(&res)),
          "HYPRE_IJVectorGetObject");
        return res;
      }

      /**
       * @brief Copies a host vector to the memory of the backend.
       * @param host Host vector.
       * @param device Array of the same size.
       */
      static void upload(const Math::Vector<Real>& host, const Array<Real>& device)
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
      static void download(const Array<Real>& device, Math::Vector<Real>& host)
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
  };
}

#endif
