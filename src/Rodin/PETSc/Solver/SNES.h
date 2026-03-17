#ifndef RODIN_SOLVER_PETSC_SNES_H
#define RODIN_SOLVER_PETSC_SNES_H

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
  class SNES
    : public PETSc::Object<::SNES>, public NewtonSolverBase<KSP>
  {
    public:
      using HandleType = ::SNES;

      using LinearSystemType = PETSc::Math::LinearSystem;

      using VectorType = ::Vec;

      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      using PetscParent = PETSc::Object<HandleType>;

      using NewtonSolverParent = NewtonSolverBase<KSP>;

      using NewtonSolverParent::solve;

      /**
       * @brief Construct SNES from a Rodin KSP linear solver.
       * @param ksp Rodin KSP wrapper; the associated ProblemBase is
       *   obtained via ksp.getProblem().
       */
      explicit SNES(KSP& ksp);

      virtual ~SNES() override;

      SNES& setType(::SNESType type) noexcept;

      SNES& setTolerances(PetscReal abstol,
                          PetscReal rtol,
                          PetscReal stol,
                          PetscInt maxIt,
                          PetscInt maxF) noexcept;

      void solve(VectorType& x) override;

      HandleType& getHandle() noexcept override;

      const HandleType& getHandle() const noexcept override;

      virtual SNES* copy() const noexcept override
      {
        return new SNES(*this);
      }

    private:
      static PetscErrorCode Residual(::SNES snes, ::Vec x, ::Vec f, void* ctx);
      static PetscErrorCode Jacobian(::SNES snes, ::Vec x, ::Mat J, ::Mat P, void* ctx);

    private:
      HandleType m_snes;
      ::SNESType m_type;
      PetscReal m_abstol, m_rtol, m_stol;
      PetscInt m_maxIt, m_maxF;
  };
}

namespace Rodin::PETSc::Solver
{
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_SOLVER_PETSC_SNES_H
