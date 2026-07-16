/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FDProbe.h
 * @brief Finite-difference probes for assembled Newton systems.
 */
#ifndef RODIN_TEST_FDPROBE_H
#define RODIN_TEST_FDPROBE_H

#include <algorithm>
#include <cmath>
#include <utility>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/Problem.h"

namespace Rodin::Test
{
  /**
   * @brief Report produced by a finite-difference tangent probe.
   */
  struct FDProbeReport
  {
    Real epsilon = 0.0;
    Real absoluteError = 0.0;
    Real relativeError = 0.0;
    Real tangentNorm = 0.0;
    Real finiteDifferenceNorm = 0.0;
  };

  /**
   * @brief Finite-difference probe for a nonlinear assembled problem.
   *
   * The probe checks the Newton consistency identity
   * @f[
   *   A(u)w \approx
   *   -\frac{b(u + \epsilon w) - b(u - \epsilon w)}{2\epsilon},
   * @f]
   * where Rodin's @c Problem stores the Newton system as
   * @f$ A(u)\Delta u = b(u) = -R(u) @f$.
   *
   * For problems with essential Dirichlet conditions, use the overload taking
   * an explicit direction and pass an admissible direction, i.e. one whose
   * constrained DOFs are zero.
   *
   * @tparam ProblemType Rodin variational problem type.
   */
  template <class ProblemType>
  class FDProbe
  {
    public:
      using Problem = ProblemType;
      using LinearSystemType =
        typename Problem::LinearSystemType;
      using VectorType =
        typename LinearSystemType::VectorType;

      explicit FDProbe(Problem& problem)
        : m_problem(problem)
      {}

      /**
       * @brief Tests the problem with a deterministic unit perturbation.
       * @param epsilon Central-difference step.
       * @return Error report.
       */
      FDProbeReport test(Real epsilon = 1e-6)
      {
        return test(makeDeterministicDirection(), epsilon);
      }

      /**
       * @brief Tests the problem in a caller-supplied perturbation direction.
       * @param direction Direction @f$w@f$ used in the finite difference.
       * @param epsilon Central-difference step.
       * @return Error report.
       */
      FDProbeReport test(const VectorType& direction, Real epsilon = 1e-6)
      {
        auto& state = m_problem.getTrialFunction().getSolution();
        const auto stateData = state.getData();

        state.getData() = stateData + epsilon * direction;
        m_problem.assemble();
        const auto rhsPlus = m_problem.getLinearSystem().getVector();

        state.getData() = stateData - epsilon * direction;
        m_problem.assemble();
        const auto rhsMinus = m_problem.getLinearSystem().getVector();

        state.getData() = stateData;
        m_problem.assemble();
        const auto tangent = m_problem.getLinearSystem().getOperator();

        const auto finiteDifference =
          -(1.0 / (2.0 * epsilon)) * (rhsPlus - rhsMinus);
        const auto tangentAction = tangent * direction;
        const auto error = tangentAction - finiteDifference;

        FDProbeReport report;
        report.epsilon = epsilon;
        report.absoluteError = error.norm();
        report.tangentNorm = tangentAction.norm();
        report.finiteDifferenceNorm = finiteDifference.norm();
        report.relativeError =
          report.absoluteError
          / std::max<Real>(
              1.0,
              std::max(report.tangentNorm, report.finiteDifferenceNorm));
        return report;
      }

    private:
      VectorType makeDeterministicDirection()
      {
        const auto& state = m_problem.getTrialFunction().getSolution();
        auto direction = state.getData();
        for (decltype(direction.size()) i = 0; i < direction.size(); ++i)
        {
          const Real x = static_cast<Real>(i + 1);
          direction[i] =
            0.5 * std::sin(0.37 * x)
            + 0.25 * std::cos(0.19 * x);
        }

        const Real norm = direction.norm();
        if (norm > 0.0)
          direction *= 1.0 / norm;
        return direction;
      }

      Problem& m_problem;
  };

  /**
   * @brief Finite-difference probe for an assembled linear system.
   *
   * This specialization probes the affine residual
   * @f[
   *   F(x) = Ax - b
   * @f]
   * and checks
   * @f[
   *   Aw \approx \frac{F(x + \epsilon w) - F(x - \epsilon w)}{2\epsilon}.
   * @f]
   */
  template <class Operator, class Vector>
  class FDProbe<Math::LinearSystem<Operator, Vector>>
  {
    public:
      using LinearSystemType = Math::LinearSystem<Operator, Vector>;
      using VectorType = Vector;

      explicit FDProbe(LinearSystemType& system)
        : m_system(system)
      {}

      FDProbeReport test(Real epsilon = 1e-6)
      {
        return test(makeDeterministicDirection(), epsilon);
      }

      FDProbeReport test(const VectorType& direction, Real epsilon = 1e-6)
      {
        auto& solution = m_system.getSolution();
        if (solution.size() != direction.size())
        {
          solution.resize(direction.size());
          solution.setZero();
        }

        const auto solutionData = solution;

        solution = solutionData + epsilon * direction;
        VectorType residualPlus =
          m_system.getOperator() * solution - m_system.getVector();

        solution = solutionData - epsilon * direction;
        VectorType residualMinus =
          m_system.getOperator() * solution - m_system.getVector();

        solution = solutionData;

        const auto finiteDifference =
          (1.0 / (2.0 * epsilon)) * (residualPlus - residualMinus);
        const auto tangentAction = m_system.getOperator() * direction;
        const auto error = tangentAction - finiteDifference;

        FDProbeReport report;
        report.epsilon = epsilon;
        report.absoluteError = error.norm();
        report.tangentNorm = tangentAction.norm();
        report.finiteDifferenceNorm = finiteDifference.norm();
        report.relativeError =
          report.absoluteError
          / std::max<Real>(
              1.0,
              std::max(report.tangentNorm, report.finiteDifferenceNorm));
        return report;
      }

    private:
      VectorType makeDeterministicDirection()
      {
        VectorType direction(m_system.getOperator().cols());
        for (decltype(direction.size()) i = 0; i < direction.size(); ++i)
        {
          const Real x = static_cast<Real>(i + 1);
          direction[i] =
            0.5 * std::sin(0.37 * x)
            + 0.25 * std::cos(0.19 * x);
        }

        const Real norm = direction.norm();
        if (norm > 0.0)
          direction *= 1.0 / norm;
        return direction;
      }

      LinearSystemType& m_system;
  };
}

#endif
