#ifndef RODIN_SOLVER_PETSC_SNES_H
#define RODIN_SOLVER_PETSC_SNES_H

#include <petscsnes.h>

#include "Rodin/PETSc/Object.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "KSP.h"

namespace Rodin::Solver
{
  class SNES
    : public PETSc::Object<::SNES>, public NewtonSolverBase<PETSc::Math::LinearSystem, ::KSP>
  {
    public:
      using HandleType = ::SNES;
      using KSPType = ::KSP;
      using MatrixType = ::Mat;
      using VectorType = ::Vec;
      using LinearSystemType = PETSc::Math::LinearSystem;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using PetscParent = PETSc::Object<HandleType>;
      using NewtonParent = NewtonSolverBase<LinearSystemType, KSPType>;
      using NewtonParent::solve;

      /**
       * @brief Construct SNES from a variational problem.
       * @param pb Variational problem associated to this Newton solver.
       * @param comm PETSc communicator used to create the SNES handle.
       */
      explicit SNES(ProblemBaseType& pb, MPI_Comm comm = PETSC_COMM_WORLD);

      virtual ~SNES() override;

      SNES& setType(::SNESType type) noexcept;

      SNES& setTolerances(PetscReal abstol,
                          PetscReal rtol,
                          PetscReal stol,
                          PetscInt maxIt,
                          PetscInt maxF) noexcept;

      SNES& setKSP(KSPType ksp) noexcept;

      template <class KSPSubclass>
      SNES& setKSP(KSPSubclass& ksp) noexcept
      {
        return setKSP(ksp.getHandle());
      }

      void solve(VectorType b, VectorType x);
      void solve(VectorType& x) override;
      void solve(LinearSystemType& system);

      const KSPType& getLinearSolver() const override;

      KSPType& getLinearSolver() override;

      HandleType& getHandle() noexcept override;

      const HandleType& getHandle() const noexcept override;

      virtual SNES* copy() const noexcept override
      {
        return new SNES(*this);
      }

    private:
      static PetscErrorCode assembleResidual(::SNES snes, ::Vec x, ::Vec f, void* ctx);
      static PetscErrorCode assembleJacobian(::SNES snes, ::Vec x, ::Mat J, ::Mat P, void* ctx);

    private:
      HandleType m_snes;
      ::SNESType m_type;
      PetscReal m_abstol, m_rtol, m_stol;
      PetscInt m_maxIt, m_maxF;
      KSPType m_kspHandle;
  };
}

namespace Rodin::PETSc::Solver
{
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_SOLVER_PETSC_SNES_H
