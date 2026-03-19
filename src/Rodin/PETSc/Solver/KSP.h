/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_PETSC_KSP_H
#define RODIN_SOLVER_PETSC_KSP_H

/**
 * @file
 * @brief PETSc KSP linear solver wrapper for Rodin variational problems.
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
   * @brief PETSc KSP (Krylov) linear solver wrapper.
   *
   * Inherits LinearSolverBase<Mat,Vec,PetscScalar> for the generic interface,
   * and PETSc::Object for automatic cleanup of any forgotten handles.
   *
   * Combines programmatic configuration with command‐line overrides.
   */
  class KSP
    : public LinearSolverBase<PETSc::Math::LinearSystem>, public PETSc::Object<::KSP>
  {
    public:
      using HandleType = ::KSP;
      using ScalarType   = PetscScalar;
      using OperatorType = ::Mat;
      using VectorType   = ::Vec;
      using LinearSystemType = PETSc::Math::LinearSystem;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
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
