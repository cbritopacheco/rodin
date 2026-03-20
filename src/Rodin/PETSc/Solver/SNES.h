#ifndef RODIN_SOLVER_PETSC_SNES_H
#define RODIN_SOLVER_PETSC_SNES_H

/**
 * @file SNES.h
 * @brief PETSc SNES (Scalable Nonlinear Equations Solvers) wrapper for Rodin.
 *
 * Wraps the PETSc `SNES` context for solving nonlinear systems
 * @f$ F(x) = 0 @f$ using Newton-type methods.  At each Newton iteration
 * the Jacobian system is solved by the associated
 * @ref Rodin::Solver::KSP linear solver.
 *
 * @see Rodin::Solver::KSP,
 *      Rodin::Solver::NewtonSolverBase,
 *      Rodin::PETSc::Math::LinearSystem
 */

#include <petscsnes.h>

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
    using LinearSystemType = PETSc::Math::LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @brief PETSc SNES (Scalable Nonlinear Equations Solvers) wrapper.
   *
   * Wraps the PETSc `SNES` context for solving nonlinear systems
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
      /// @brief Handle type for the raw PETSc `SNES` context pointer.
      using HandleType = ::SNES;

      /// @brief Linear system type coupling @f$ A @f$, @f$ \mathbf{b} @f$, and @f$ \mathbf{x} @f$.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief PETSc vector type (`::Vec`) for the nonlinear residual and solution.
      using VectorType = ::Vec;

      /// @brief Base problem type that provides the linear system.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// @brief Parent class providing PETSc object handle management.
      using PetscParent = PETSc::Object<HandleType>;

      /// @brief Parent class providing the Newton solver interface.
      using NewtonSolverParent = NewtonSolverBase<KSP>;

      using NewtonSolverParent::solve;

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
       * @brief Solves the nonlinear system @f$ F(x) = 0 @f$.
       * @param[in,out] x Initial guess on input; solution on output.
       */
      void solve(VectorType& x) override;

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
      static PetscErrorCode Residual(::SNES snes, ::Vec x, ::Vec f, void* ctx);
      static PetscErrorCode Jacobian(::SNES snes, ::Vec x, ::Mat J, ::Mat P, void* ctx);

    private:
      HandleType m_snes;  ///< Underlying PETSc SNES context.
      ::SNESType m_type;  ///< Requested SNES algorithm type.
      PetscReal m_abstol, ///< Absolute convergence tolerance.
                m_rtol,   ///< Relative convergence tolerance.
                m_stol;   ///< Step norm convergence tolerance.
      PetscInt m_maxIt,   ///< Maximum nonlinear iterations.
               m_maxF;    ///< Maximum function evaluations.
  };
}

namespace Rodin::PETSc::Solver
{
  /**
   * @brief PETSc namespace alias to @ref Rodin::Solver::SNES.
   */
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_SOLVER_PETSC_SNES_H
