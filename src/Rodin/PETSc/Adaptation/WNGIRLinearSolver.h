/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ADAPTATION_WNGIRLINEARSOLVER_H
#define RODIN_PETSC_ADAPTATION_WNGIRLINEARSOLVER_H

#include <algorithm>
#include <cstring>
#include <petscksp.h>

#include "Rodin/Adaptation/WNGIRParameters.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/PETSc/Solver/KSP.h"
#include "Rodin/Types.h"

namespace Rodin::PETSc::Adaptation
{
  template <class ProblemT>
  class WNGIRLinearSolver
  {
    public:
      explicit WNGIRLinearSolver(ProblemT& problem)
        : m_ksp(problem)
      {
        m_ksp.setPrefix("wngir_");
      }

      bool solve(
          ProblemT& problem,
          Rodin::Math::Vector<Real>& solution,
          std::size_t& iterations,
          Real& error,
          const Rodin::Adaptation::WNGIRParameters& parameters)
      {
        auto& system = problem.getLinearSystem();
        const ::Mat A = system.getOperator();
        PetscObjectState state = 0;
        if (MatGetNonzeroState(A, &state) != PETSC_SUCCESS)
          return false;

        const bool patternChanged = m_hasSolved && state != m_nonzeroState;
        if (m_hasSolved)
          setFactorReuse(!patternChanged);

        m_ksp.setTolerances(
            parameters.cgRelativeTolerance > Real(0)
              ? parameters.cgRelativeTolerance : Real(1e-6),
            PETSC_DEFAULT, PETSC_DEFAULT,
            parameters.cgMaxIterations > 0
              ? static_cast<PetscInt>(parameters.cgMaxIterations) : 2000);
        problem.solve(m_ksp);

        ::KSPConvergedReason reason;
        PetscInt its = 0;
        PetscReal residual = 0, rhsNorm = 0;
        KSPGetConvergedReason(m_ksp.getHandle(), &reason);
        KSPGetIterationNumber(m_ksp.getHandle(), &its);
        KSPGetResidualNorm(m_ksp.getHandle(), &residual);
        VecNorm(system.getVector(), NORM_2, &rhsNorm);
        iterations = static_cast<std::size_t>(std::max<PetscInt>(its, 0));
        error = residual / std::max<Real>(rhsNorm, Real(1e-30));

        const ::Vec x = system.getSolution();
        PetscInt n = 0;
        VecGetSize(x, &n);
        solution.resize(n);
        const PetscScalar* values = nullptr;
        VecGetArrayRead(x, &values);
        for (PetscInt i = 0; i < n; ++i)
          solution(i) = PetscRealPart(values[i]);
        VecRestoreArrayRead(x, &values);

        m_nonzeroState = state;
        m_hasSolved = true;
        setFactorReuse(true);
        return reason > 0 && solution.allFinite();
      }

      Rodin::Solver::KSP& getKSP() noexcept
      {
        return m_ksp;
      }

    private:
      void setFactorReuse(bool reuse)
      {
        ::PC pc = nullptr;
        KSPGetPC(m_ksp.getHandle(), &pc);
        const char* type = nullptr;
        PCGetType(pc, &type);
        if (type && (!std::strcmp(type, PCLU)
                  || !std::strcmp(type, PCCHOLESKY)))
        {
          PCFactorSetReuseOrdering(pc, reuse ? PETSC_TRUE : PETSC_FALSE);
          PCFactorSetReuseFill(pc, reuse ? PETSC_TRUE : PETSC_FALSE);
        }
      }

      Rodin::Solver::KSP m_ksp;
      PetscObjectState m_nonzeroState = 0;
      bool m_hasSolved = false;
  };
}

#endif
