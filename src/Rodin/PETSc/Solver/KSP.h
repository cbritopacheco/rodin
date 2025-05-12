#ifndef RODIN_SOLVER_PETSC_KSP_H
#define RODIN_SOLVER_PETSC_KSP_H

#include <petscksp.h>

#include "Rodin/Solver/CG.h"

namespace Rodin::Solver
{
  class KSP : public SolverBase<::Mat, ::Vec, PetscScalar>
  {
    public:
      using ScalarType = PetscScalar;
      using VectorType = ::Vec;
      using OperatorType = ::Mat;
      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;
      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;
      using Parent::solve;

      KSP(ProblemType& pb)
        : Parent(pb)
      {}

      KSP(const KSP& other)
        : Parent(other)
      {}

      KSP(KSP&& other)
        : Parent(std::move(other))
      {}

      ~KSP() = default;

      virtual void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        KSPSetType(m_ksp, m_type);
        KSPSetTolerances(m_ksp, m_rtol, m_abstol, m_dtol, m_maxIt);
        if (m_preconditioner)
          KSPSetOperators(m_ksp, A, m_preconditioner.value().get());
        else
          KSPSetOperators(m_ksp, A, A);
        KSPSolve(m_ksp, b, x);
      }

      KSP& create()
      {
        KSPCreate(PETSC_COMM_WORLD, &m_ksp);
        return *this;
      }

      KSP& setType(::KSPType type)
      {
        m_type = type;
        return *this;
      }

      KSP& setTolerances(PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxIt)
      {
        m_rtol = rtol;
        m_abstol = abstol;
        m_dtol = dtol;
        m_maxIt = maxIt;
        return *this;
      }

      KSP& setPreconditioner(::Mat& preconditioner)
      {
        m_preconditioner = preconditioner;
        return *this;
      }

      KSP& destroy()
      {
        KSPDestroy(&m_ksp);
        return *this;
      }

      virtual KSP* copy() const noexcept override
      {
        return new KSP(*this);
      }

    private:
      ::KSP m_ksp;
      ::KSPType m_type;
      PetscReal m_rtol;
      PetscReal m_abstol;
      PetscReal m_dtol;
      PetscInt m_maxIt;
      std::optional<std::reference_wrapper<::Mat>> m_preconditioner;
  };
}

#endif

