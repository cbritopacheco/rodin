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
      // (snes, x, f, ctx): assemble residual f = F(x)
      using FunctionCallbackType = PetscErrorCode (*)(::SNES, ::Vec, ::Vec, void*);
      // (snes, x, J, P, ctx): assemble Jacobian/preconditioner operators
      using JacobianCallbackType = PetscErrorCode (*)(::SNES, ::Vec, ::Mat, ::Mat, void*);
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

      SNES& setFunction(
          FunctionCallbackType f,
          void* ctx = PETSC_NULLPTR,
          VectorType residual = PETSC_NULLPTR);

      SNES& setFunction(
          FunctionCallbackType f,
          VectorType residual);

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

      SNES& setJacobian(JacobianCallbackType j,
                        MatrixType jacobian,
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
