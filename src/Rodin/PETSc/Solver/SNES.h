/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_SOLVER_SNES_H
#define RODIN_PETSC_SOLVER_SNES_H

/**
 * @file SNES.h
 * @brief PETSc SNES (Scalable Nonlinear Equations Solvers) wrapper for Rodin.
 *
 * Wraps the PETSc @c SNES context for solving nonlinear systems
 * @f$ F(x) = 0 @f$ using Newton-type methods.  At each Newton iteration
 * the Jacobian system is solved by the associated
 * @ref Rodin::Solver::KSP linear solver.
 *
 * @see Rodin::Solver::KSP,
 *      Rodin::Solver::NewtonSolverBase,
 *      Rodin::PETSc::Math::LinearSystem
 */

#include <petscsnes.h>

#include <functional>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/PETSc/Object.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "KSP.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Traits specialization for Solver::KSP (Rodin wrapper for PETSc KSP).
   *
   * Maps Solver::KSP to the PETSc linear system type so that
   * NewtonSolverBase<Solver::KSP> can deduce its LinearSystemType.
   */
  template <>
  struct Traits<Solver::KSP>
  {
    /// @brief Linear system type.
      using LinearSystemType = PETSc::Math::LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @brief PETSc SNES (Scalable Nonlinear Equations Solvers) wrapper.
   *
   * Wraps the PETSc @c SNES context for solving nonlinear systems
   * @f$ F(x) = 0 @f$ using Newton-type methods.  Inherits
   * @ref Rodin::Solver::NewtonSolverBase<KSP> so that the linear
   * sub-problems arising at each Newton step are solved by the
   * associated @ref Rodin::Solver::KSP solver, and
   * @ref Rodin::PETSc::Object<::SNES> for automatic handle cleanup.
   *
   * Supports both programmatic configuration (`setType`, `setTolerances`)
   * and PETSc command-line overrides (`-snes_type`, `-snes_rtol`, …).
   *
   * @see Rodin::Solver::KSP,
   *      Rodin::Solver::NewtonSolverBase,
   *      Rodin::PETSc::Math::LinearSystem
   */
  class SNES
    : public PETSc::Object<::SNES>, public NewtonSolverBase<KSP>
  {
    public:
      /// @brief Handle type for the raw PETSc @c SNES context pointer.
      using HandleType = ::SNES;

      /// @brief Linear system type coupling @f$ A @f$, @f$ \mathbf{b} @f$, and @f$ \mathbf{x} @f$.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief PETSc vector type (@c Vec) for the nonlinear residual and solution.
      using VectorType = ::Vec;

      /**
       * @brief Callback invoked before residual/Jacobian assembly to synchronize
       *        state-dependent Rodin fields from the current SNES iterate.
       *
       * The argument is the current iterate expressed as a Rodin
       * @ref Rodin::PETSc::Math::Vector (a PETSc @c Vec).  Typical usage:
       * @code
       * snes.setStateUpdate([&](const PETSc::Math::Vector& x) {
       *   uState.setData(x, 0);
       *   pState.setData(x, uh.getSize());
       * });
       * @endcode
       */
      using StateUpdate = std::function<void(const PETSc::Math::Vector&)>;

      /// @brief Base problem type that provides the linear system.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// @brief Parent class providing PETSc object handle management.
      using PetscParent = PETSc::Object<HandleType>;

      /// @brief Parent class providing the Newton solver interface.
      using NewtonSolverParent = NewtonSolverBase<KSP>;

      using NewtonSolverParent::solve;

      /// @brief Default nonlinear relative residual tolerance.
      static constexpr PetscReal DEFAULT_RTOL   = 1e-8;
      /// @brief Default nonlinear absolute residual tolerance.
      static constexpr PetscReal DEFAULT_ABSTOL = 1e-50;
      /// @brief Default step norm convergence tolerance.
      static constexpr PetscReal DEFAULT_STOL   = 1e-8;
      /// @brief Default maximum number of nonlinear iterations.
      static constexpr PetscInt  DEFAULT_MAXIT  = 50;
      /// @brief Default maximum number of residual evaluations.
      static constexpr PetscInt  DEFAULT_MAXF   = 10000;

      /**
       * @brief Construct SNES from a Rodin KSP linear solver.
       * @param ksp Rodin KSP wrapper; the associated ProblemBase is
       *   obtained via ksp.getProblem().
       */
      explicit SNES(KSP& ksp);

      virtual ~SNES() override;

      /**
       * @brief Selects the nonlinear solver algorithm.
       * @param type PETSc SNES type string (e.g. `SNESNEWTONLS`).
       * @returns Reference to `*this`.
       */
      SNES& setType(::SNESType type) noexcept;

      /**
       * @brief Sets nonlinear convergence tolerances and limits.
       * @param abstol Absolute convergence tolerance.
       * @param rtol   Relative convergence tolerance.
       * @param stol   Convergence tolerance in terms of step norm.
       * @param maxIt  Maximum number of nonlinear iterations.
       * @param maxF   Maximum number of function evaluations.
       * @returns Reference to `*this`.
       */
      SNES& setTolerances(PetscReal abstol,
                          PetscReal rtol,
                          PetscReal stol,
                          PetscInt maxIt,
                          PetscInt maxF) noexcept;

      /**
       * @brief Sets an optional callback invoked before residual/Jacobian assembly.
       *
       * The callback receives the current SNES iterate as a Rodin
       * @ref PETSc::Math::Vector.  It is responsible for copying the relevant
       * sub-vectors back into the GridFunctions that appear in the nonlinear
       * variational form.
       *
       * @code
       * snes.setStateUpdate([&](const PETSc::Math::Vector& x) {
       *   uState.setData(x, 0);
       *   pState.setData(x, uh.getSize());
       * });
       * @endcode
       *
       * @param update Callback to invoke on each residual/Jacobian assembly.
       * @returns Reference to `*this`.
       */
      SNES& setStateUpdate(StateUpdate update);

      /**
       * @brief Solves the nonlinear system @f$ F(x) = 0 @f$.
       * @param[in,out] x Initial guess on input; solution on output.
       */
      void solve(VectorType& x) override;

      /**
       * @brief Solves the nonlinear system using the linear system's solution
       *        vector as both the initial guess and the solution storage.
       *
       * The solution vector is obtained directly from
       * @c ksp.getProblem().getLinearSystem().getSolution(), so no manual
       * initial-guess packing is needed. A full @c Problem::assemble() seeds
       * this vector from the trial function data, so the initial iterate is
       * the trial functions' state — consistent with
       * @c NewtonSolverBase::solve(GridFunction&). After the first solve the
       * solution vector retains its value, providing a natural warm start for
       * subsequent time steps.
       */
      void solve();

      /**
       * @brief Returns the number of SNES iterations from the most recent solve.
       * @returns PETSc iteration count.
       */
      PetscInt getIterationNumber() const;

      /**
       * @brief Returns `true` if the most recent @ref solve() converged.
       *
       * Wraps @c SNESGetConvergedReason: a positive reason indicates
       * convergence, a negative reason indicates divergence.
       *
       * @returns `true` when converged, `false` when diverged or not yet solved.
       */
      bool converged() const;

      /**
       * @brief Returns PETSc's convergence reason from the most recent solve.
       * @returns Positive values for convergence, negative values for
       *          divergence, and zero before PETSc has decided.
       */
      ::SNESConvergedReason getConvergedReason() const;

      /// @brief Returns a mutable reference to the underlying PETSc SNES handle.
      /// @returns Mutable reference to the SNES handle.
      HandleType& getHandle() noexcept override;

      /// @brief Returns a read-only reference to the underlying PETSc SNES handle.
      /// @returns Const reference to the SNES handle.
      const HandleType& getHandle() const noexcept override;

      /// @brief Creates a heap-allocated copy of this SNES solver.
      /// @returns Pointer to the cloned SNES instance.
      virtual SNES* copy() const noexcept override
      {
        return new SNES(*this);
      }

    private:
      static PetscErrorCode Update(::Vec x, void* ctx);
      static PetscErrorCode Assemble(
        ::Vec x, void* ctx, Variational::AssemblyTarget target);
      static PetscErrorCode Residual(::SNES snes, ::Vec x, ::Vec f, void* ctx);
      static PetscErrorCode Jacobian(::SNES snes, ::Vec x, ::Mat J, ::Mat P, void* ctx);

    private:
      HandleType m_snes;   ///< Underlying PETSc SNES context.
      ::SNESType m_type;   ///< Requested SNES algorithm type.
      ::PetscReal m_abstol, ///< Absolute convergence tolerance.
        m_rtol, ///< Relative convergence tolerance.
        m_stol; ///< Step norm convergence tolerance.
      ::PetscInt m_maxIt, ///< Maximum nonlinear iterations.
        m_maxF; ///< Maximum function evaluations.
      StateUpdate m_update; ///< Optional state synchronization callback.
      Optional<::PetscObjectState> m_lhsAssembled;
      Optional<::PetscObjectState> m_rhsAssembled;
      Optional<::PetscObjectState> m_updated;
  };
}

namespace Rodin::PETSc::Solver
{
  /**
   * @brief PETSc namespace alias to @ref Rodin::Solver::SNES.
   */
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_PETSC_SOLVER_SNES_H
