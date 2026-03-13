#ifndef RODIN_SOLVER_PETSC_SNES_H
#define RODIN_SOLVER_PETSC_SNES_H

#include <functional>
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
      using FunctionType =
        std::function<PetscErrorCode(HandleType, VectorType, VectorType)>;
      using JacobianType =
        // (snes, x, A, P): A is Jacobian operator, P is preconditioner operator.
        std::function<PetscErrorCode(HandleType, VectorType, MatrixType, MatrixType)>;
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

      SNES& setFunction(FunctionType f, VectorType residual = PETSC_NULLPTR);

      SNES& setJacobian(JacobianType j,
                        MatrixType jacobian = PETSC_NULLPTR,
                        MatrixType preconditioner = PETSC_NULLPTR);

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
      FunctionType m_function;
      JacobianType m_jacobian;
      VectorType m_residual;
      MatrixType m_jacobianOperator;
      MatrixType m_preconditionerOperator;

      static PetscErrorCode FunctionCallback(
          ::SNES snes, ::Vec x, ::Vec f, void* ctx);

      static PetscErrorCode JacobianCallback(
          ::SNES snes, ::Vec x, ::Mat A, ::Mat P, void* ctx);
  };
}

namespace Rodin::PETSc::Solver
{
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_SOLVER_PETSC_SNES_H
