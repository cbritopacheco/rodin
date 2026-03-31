/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_PETSC_KSP_H
#define RODIN_SOLVER_PETSC_KSP_H

/**
 * @file KSP.h
 * @brief PETSc KSP (Krylov subspace) linear solver wrapper for Rodin.
 *
 * Wraps the PETSc `KSP` context and implements the Rodin
 * @ref Rodin::Solver::LinearSolverBase interface so that any PETSc
 * Krylov solver (CG, GMRES, BiCGStab, …) can be used to solve the
 * linear system @f$ A\mathbf{x} = \mathbf{b} @f$ assembled by a
 * Rodin variational problem.
 *
 * The class supports:
 * - Programmatic configuration via `setType()`, `setTolerances()`,
 *   `setPreconditioner()`.
 * - Command-line overrides via `KSPSetFromOptions` (called inside `solve`).
 * - Automatic solution vector allocation when `x == PETSC_NULL`.
 *
 * @see Rodin::PETSc::Solver::CG,
 *      Rodin::PETSc::Solver::GMRES,
 *      Rodin::PETSc::Solver::SNES
 */

#include <petscksp.h>
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/PETSc/Object.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @brief PETSc KSP (Krylov subspace) linear solver wrapper.
   *
   * Wraps the PETSc `KSP` context and inherits both
   * @ref Rodin::Solver::LinearSolverBase (for the generic solver interface)
   * and @ref Rodin::PETSc::Object (for automatic handle cleanup).
   *
   * Combines programmatic configuration (tolerances, type, preconditioner)
   * with PETSc command-line overrides (`-ksp_type`, `-ksp_rtol`, …).
   *
   * @see Rodin::Solver::CG<PETSc::Math::LinearSystem>,
   *      Rodin::Solver::GMRES<PETSc::Math::LinearSystem>,
   *      Rodin::Solver::SNES
   */
  class KSP
    : public LinearSolverBase<PETSc::Math::LinearSystem>, public PETSc::Object<::KSP>
  {
    public:
      /// @brief Handle type for the raw PETSc `KSP` context pointer.
      using HandleType = ::KSP;
      /// @brief Scalar type (`PetscScalar`) used for residual norms and tolerances.
      using ScalarType   = PetscScalar;
      /// @brief PETSc matrix type (`::Mat`) for the system operator and preconditioner.
      using OperatorType = ::Mat;
      /// @brief PETSc vector type (`::Vec`) for the right-hand side and solution.
      using VectorType   = ::Vec;
      /// @brief Linear system type coupling @f$ A @f$, @f$ \mathbf{b} @f$, and @f$ \mathbf{x} @f$.
      using LinearSystemType = PETSc::Math::LinearSystem;
      /// @brief Base problem type that provides the linear system to solve.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Parent class providing the generic `LinearSolverBase` interface.
      using Parent = LinearSolverBase<LinearSystemType>;
      using Parent::solve;
      /**
       * @brief Construct and create the PETSc KSP object.
       *
       * Initializes to PETSC defaults.
       *
       * @param pb   Variational problem this solver will solve.
       */
      explicit KSP(ProblemBaseType& pb);

      virtual ~KSP() override;

      /**
       * @brief Solve @f$ Ax = b @f$, allocating @f$ x @f$ if null.
       *
       * If `x == PETSC_NULL`, automatically `VecDuplicate(b, &x)` and zero it.
       * Otherwise uses `x` contents as initial guess.
       *
       * Applies programmatic settings, then SetFromOptions, then KSPSolve.
       *
       * @param b Linear system containing matrix and right-hand side vector.
       */
      void solve(LinearSystemType& b) override;

      /**
       * @brief Selects the Krylov subspace method.
       * @param type PETSc KSP type string (e.g. `KSPCG`, `KSPGMRES`).
       * @returns Reference to `*this` for method chaining.
       */
      KSP& setType(::KSPType type) noexcept;

      /**
       * @brief Sets convergence tolerances and maximum iteration count.
       * @param rtol   Relative decrease in residual norm.
       * @param abstol Absolute residual norm threshold.
       * @param dtol   Divergence tolerance.
       * @param maxIt  Maximum number of iterations.
       * @returns Reference to `*this`.
       */
      KSP& setTolerances(PetscReal rtol,
                         PetscReal abstol,
                         PetscReal dtol,
                         PetscInt  maxIt) noexcept;

      /**
       * @brief Sets an explicit preconditioner matrix.
       * @param P PETSc matrix to use as preconditioner operator.
       * @returns Reference to `*this`.
       */
      KSP& setPreconditioner(OperatorType P) noexcept;

      KSP& setPrefix(const Optional<std::string>& prefix) noexcept;

      /// @brief Returns a mutable reference to the underlying PETSc KSP handle.
      /// @returns Mutable reference to the KSP handle.
      HandleType& getHandle() noexcept override;

      /// @brief Returns a read-only reference to the underlying PETSc KSP handle.
      /// @returns Const reference to the KSP handle.
      const HandleType& getHandle() const noexcept override;

      /// @brief Creates a heap-allocated copy of this solver.
      /// @returns Pointer to the cloned KSP instance.
      virtual KSP* copy() const noexcept override
      {
        return new KSP(*this);
      }

    private:
      HandleType   m_ksp;  ///< Underlying PETSc KSP context.
      ::KSPType    m_type; ///< Requested KSP algorithm type.
      PetscReal    m_rtol, ///< Relative convergence tolerance.
                   m_abstol, ///< Absolute convergence tolerance.
                   m_dtol; ///< Divergence tolerance.
      PetscInt     m_maxIt; ///< Maximum iteration count.
      std::optional<OperatorType> m_preconditioner; ///< Optional preconditioner matrix.
      Optional<std::string> m_prefix; ///< Optional prefix for PETSc options.
  };
} // namespace Rodin::Solver

namespace Rodin::PETSc::Solver
{
  /**
   * @brief PETSc namespace alias to @ref Rodin::Solver::KSP.
   */
  using KSP = Rodin::Solver::KSP;
}

#endif // RODIN_SOLVER_PETSC_KSP_H
