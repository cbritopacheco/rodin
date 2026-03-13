#ifndef RODIN_SOLVER_PETSC_SNES_H
#define RODIN_SOLVER_PETSC_SNES_H

#include <petscsnes.h>

#include "Rodin/PETSc/Object.h"
#include "KSP.h"

namespace Rodin::Solver
{
  class SNES : public PETSc::Object<::SNES>
  {
    public:
      using HandleType = ::SNES;
      using KSPType = ::KSP;
      using MatrixType = ::Mat;
      using VectorType = ::Vec;
      using FunctionCallbackType = ::SNESFunctionFn*;
      using JacobianCallbackType = ::SNESJacobianFn*;
      using Parent = PETSc::Object<HandleType>;

      explicit SNES(MPI_Comm comm = PETSC_COMM_WORLD);

      virtual ~SNES() override;

      SNES& setType(::SNESType type) noexcept;

      SNES& setTolerances(PetscReal abstol,
                          PetscReal rtol,
                          PetscReal stol,
                          PetscInt maxIt,
                          PetscInt maxF) noexcept;

      SNES& setKSP(KSPType ksp) noexcept;

      SNES& setFunction(
          FunctionCallbackType f,
          void* ctx = PETSC_NULLPTR,
          VectorType residual = PETSC_NULLPTR);

      template <class Context>
      SNES& setFunction(
          FunctionCallbackType f,
          Context& ctx,
          VectorType residual = PETSC_NULLPTR)
      {
        return setFunction(f, static_cast<void*>(&ctx), residual);
      }

      SNES& setJacobian(JacobianCallbackType j,
                        void* ctx = PETSC_NULLPTR,
                        MatrixType jacobian = PETSC_NULLPTR,
                        MatrixType preconditioner = PETSC_NULLPTR);

      template <class Context>
      SNES& setJacobian(JacobianCallbackType j,
                        Context& ctx,
                        MatrixType jacobian = PETSC_NULLPTR,
                        MatrixType preconditioner = PETSC_NULLPTR)
      {
        return setJacobian(j, static_cast<void*>(&ctx), jacobian, preconditioner);
      }

      void solve(VectorType b, VectorType x);

      HandleType& getHandle() noexcept override;

      const HandleType& getHandle() const noexcept override;

      virtual SNES* copy() const noexcept
      {
        return new SNES(*this);
      }

    private:
      HandleType m_snes;
      ::SNESType m_type;
      PetscReal m_abstol, m_rtol, m_stol;
      PetscInt m_maxIt, m_maxF;
      KSPType m_kspHandle;
      FunctionCallbackType m_functionCallback;
      JacobianCallbackType m_jacobianCallback;
      void* m_functionContext;
      void* m_jacobianContext;
      VectorType m_residual;
      MatrixType m_jacobianOperator;
      MatrixType m_preconditionerOperator;
  };
}

namespace Rodin::PETSc::Solver
{
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_SOLVER_PETSC_SNES_H
