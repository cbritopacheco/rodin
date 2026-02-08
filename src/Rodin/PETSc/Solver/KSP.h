/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_PETSC_KSP_H
#define RODIN_SOLVER_PETSC_KSP_H

#include <petscksp.h>
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Solver/Solver.h"
#include "Rodin/PETSc/Object.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @brief PETSc KSP (Krylov) linear solver wrapper.
   *
   * Inherits SolverBase<Mat,Vec,PetscScalar> for the generic interface,
   * and PETSc::Object for automatic cleanup of any forgotten handles.
   *
   * Combines programmatic configuration with command‐line overrides.
   */
  class KSP
    : public SolverBase<PETSc::Math::LinearSystem>, public PETSc::Object<::KSP>
  {
    public:
      using HandleType = ::KSP;
      using ScalarType   = PetscScalar;
      using OperatorType = ::Mat;
      using VectorType   = ::Vec;
      using LinearSystemType = PETSc::Math::LinearSystem;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = SolverBase<LinearSystemType>;
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

      KSP& setType(::KSPType type) noexcept;

      KSP& setTolerances(PetscReal rtol,
                         PetscReal abstol,
                         PetscReal dtol,
                         PetscInt  maxIt) noexcept;

      KSP& setPreconditioner(OperatorType P) noexcept;

      HandleType& getHandle() noexcept override;

      const HandleType& getHandle() const noexcept override;

      virtual KSP* copy() const noexcept override
      {
        return new KSP(*this);
      }

    private:
      HandleType   m_ksp;
      ::KSPType    m_type;
      PetscReal    m_rtol,
                   m_abstol,
                   m_dtol;
      PetscInt     m_maxIt;
      std::optional<OperatorType> m_preconditioner;
  };
} // namespace Rodin::Solver

namespace Rodin::PETSc::Solver
{
  using KSP = Rodin::Solver::KSP;
}

#endif // RODIN_SOLVER_PETSC_KSP_H
