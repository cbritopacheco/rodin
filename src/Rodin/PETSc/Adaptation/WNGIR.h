#ifndef RODIN_PETSC_ADAPTATION_WNGIR_H
#define RODIN_PETSC_ADAPTATION_WNGIR_H

#include <cstring>
#include <petscksp.h>

#include "Rodin/Adaptation/WNGIR.h"
#include "Rodin/PETSc/Solver/KSP.h"
#include "Rodin/PETSc/Variational.h"

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

  template <class Displacement>
  class WNGIR
  {
      using FESType = std::decay_t<
        decltype(std::declval<Displacement&>().getFiniteElementSpace())>;
      using TrialType = std::decay_t<decltype(
        PETSc::Variational::TrialFunction(std::declval<const FESType&>()))>;
      using TestType = std::decay_t<decltype(
        PETSc::Variational::TestFunction(std::declval<const FESType&>()))>;
      using SolutionType =
        typename FormLanguage::Traits<TrialType>::SolutionType;
      using ProblemType = Rodin::Variational::Problem<
        PETSc::Math::LinearSystem, TrialType, TestType>;
      using BilinearFormType = Rodin::Variational::BilinearForm<
        SolutionType, FESType, FESType, ::Mat>;
      using LinearSolverType = WNGIRLinearSolver<ProblemType>;

    public:
      explicit WNGIR(Displacement& displacement)
        : m_displacement(&displacement),
          m_trial(displacement.getFiniteElementSpace()),
          m_test(displacement.getFiniteElementSpace()),
          m_problem(m_trial, m_test),
          m_bulkForm(m_trial, m_test),
          m_linearSolver(m_problem)
      {}

      void setParameters(const Rodin::Adaptation::WNGIRParameters& parameters)
      {
        m_parameters = parameters;
        m_bulkFormAssembled = false;
      }

      const Rodin::Adaptation::WNGIRReport& getReport() const noexcept
      {
        return m_report;
      }

      template <class Mesh, class PhiDerived, class GradDerived>
      Rodin::Adaptation::WNGIRReport solve(
          const Mesh& mesh,
          const std::vector<Rodin::Index>& interfaceFacets,
          const Rodin::Variational::RealFunctionBase<PhiDerived>& phi,
          const Rodin::Variational::VectorFunctionBase<Real, GradDerived>& grad)
      {
        m_report = Rodin::Adaptation::solveWNGIR(
            mesh, *m_displacement, interfaceFacets, phi, grad, m_parameters,
            m_trial, m_test, m_bulkForm, m_problem, m_linearSolver,
            m_bulkFormAssembled);
        return m_report;
      }

    private:
      Displacement* m_displacement;
      TrialType m_trial;
      TestType m_test;
      ProblemType m_problem;
      BilinearFormType m_bulkForm;
      LinearSolverType m_linearSolver;
      bool m_bulkFormAssembled = false;
      Rodin::Adaptation::WNGIRParameters m_parameters;
      Rodin::Adaptation::WNGIRReport m_report;
  };

  template <class Displacement>
  WNGIR(Displacement&) -> WNGIR<Displacement>;
}

#endif
