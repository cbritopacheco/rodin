/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ADAPTATION_WNGIRSOLVER_H
#define RODIN_PETSC_ADAPTATION_WNGIRSOLVER_H

#include <type_traits>
#include <vector>

#include "Rodin/Adaptation/WNGIRSolve.h"
#include "Rodin/PETSc/Adaptation/WNGIRLinearSolver.h"
#include "Rodin/PETSc/Variational.h"

namespace Rodin::PETSc::Adaptation
{
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
