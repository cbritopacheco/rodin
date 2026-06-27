/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRSOLVER_H
#define RODIN_ADAPTATION_WNGIRSOLVER_H

#include "WNGIREigenLinearSolver.h"
#include "WNGIRSolve.h"

namespace Rodin::Adaptation
{
  template <class Displacement>
  class WNGIR
  {
      using FESType = std::decay_t<
        decltype(std::declval<Displacement&>().getFiniteElementSpace())>;
      using TrialType = std::decay_t<
        decltype(Variational::TrialFunction(std::declval<const FESType&>()))>;
      using TestType = std::decay_t<
        decltype(Variational::TestFunction(std::declval<const FESType&>()))>;
      using ProblemType = std::decay_t<decltype(Variational::Problem(
          std::declval<TrialType&>(), std::declval<TestType&>()))>;
      using BilinearFormType = std::decay_t<decltype(Variational::BilinearForm(
          std::declval<TrialType&>(), std::declval<TestType&>()))>;
      using LinearSolverType = WNGIREigenLinearSolver<ProblemType>;

    public:
      explicit WNGIR(Displacement& u)
        : m_u(&u),
          m_duStep(u.getFiniteElementSpace()),
          m_vStep(u.getFiniteElementSpace()),
          m_stepProblem(m_duStep, m_vStep),
          m_bulkForm(m_duStep, m_vStep)
      {}

      void setParameters(const WNGIRParameters& parameters)
      {
        m_parameters = parameters;
        m_bulkFormAssembled = false;
      }

      const WNGIRParameters& getParameters() const
      {
        return m_parameters;
      }

      const WNGIRReport& getReport() const
      {
        return m_report;
      }

      template <class Mesh, class PhiDerived, class GradDerived>
      WNGIRReport solve(
          const Mesh& mesh,
          const std::vector<Rodin::Index>& interfaceFacets,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad)
      {
        m_report = solveWNGIR(
            mesh, *m_u, interfaceFacets, phi, grad, m_parameters,
            m_duStep, m_vStep, m_bulkForm, m_stepProblem,
            m_linearSolver, m_bulkFormAssembled);
        return m_report;
      }

    private:
      Displacement* m_u;
      TrialType m_duStep;
      TestType m_vStep;
      ProblemType m_stepProblem;
      BilinearFormType m_bulkForm;
      LinearSolverType m_linearSolver;
      bool m_bulkFormAssembled = false;
      WNGIRParameters m_parameters;
      WNGIRReport m_report;
  };

  template <class Displacement>
  WNGIR(Displacement&) -> WNGIR<Displacement>;
}

#endif
