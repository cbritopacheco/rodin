/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TMOP_OPTIMIZER_H
#define RODIN_ADAPTATION_TMOP_OPTIMIZER_H

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "Objective.h"

namespace Rodin::Adaptation::TMOP
{
  struct OptimizerOptions
  {
    Index maxIterations = 50;
    Real gradientTolerance = 1e-8;
    Real finiteDifferenceStep = 1e-6;
    Real initialStepSize = 0.2;
    Real minStepSize = 1e-12;
  };

  struct OptimizationReport
  {
    Index iterations = 0;
    Real initialObjective = 0;
    Real finalObjective = 0;
    Real minJacobian = 0;
    Index invalidElements = 0;
    bool converged = false;
    std::string reason;
  };

  class Optimizer
  {
    public:
      /**
       * @brief Damped finite-difference optimizer for fixed-topology geometry.
       *
       * This optimizer is deliberately conservative and fixed-topology:
       * it respects GeometryNode::fixed, rejects steps with invalid sampled
       * Jacobians, and uses backtracking. It provides a fixed-topology
       * execution path for tests and small examples while the residual and
       * tangent terms provide the LinearProblem/Newton assembly path.
       */
      Optimizer(HighOrderTriangleGeometry& geometry, const Objective& objective)
        : m_geometry(geometry),
          m_objective(objective)
      {}

      Optimizer& setOptions(const OptimizerOptions& options)
      {
        m_options = options;
        return *this;
      }

      OptimizationReport optimize() const
      {
        OptimizationReport report;
        report.initialObjective = m_objective.value(m_geometry);

        if (m_objective.invalidElementCount(m_geometry) > 0)
        {
          report.finalObjective = report.initialObjective;
          report.minJacobian = m_objective.minJacobian(m_geometry);
          report.invalidElements = m_objective.invalidElementCount(m_geometry);
          report.reason = "invalid initial geometry";
          return report;
        }

        Real current = report.initialObjective;
        for (Index it = 0; it < m_options.maxIterations; ++it)
        {
          std::vector<std::array<Real, 2>> gradient(m_geometry.nodes.size());
          Real gradientNorm2 = 0;

          for (Index n = 0; n < m_geometry.nodes.size(); ++n)
          {
            gradient[n] = {{ 0, 0 }};
            if (m_geometry.nodes[n].fixed)
              continue;

            for (std::uint8_t d = 0; d < 2; ++d)
            {
              const Real old = m_geometry.nodes[n].x[d];
              m_geometry.nodes[n].x[d] = old + m_options.finiteDifferenceStep;
              const Real plus = m_objective.value(m_geometry);
              m_geometry.nodes[n].x[d] = old - m_options.finiteDifferenceStep;
              const Real minus = m_objective.value(m_geometry);
              m_geometry.nodes[n].x[d] = old;

              gradient[n][d] =
                (plus - minus) / (Real(2) * m_options.finiteDifferenceStep);
              gradientNorm2 += gradient[n][d] * gradient[n][d];
            }
          }

          const Real gradientNorm = std::sqrt(gradientNorm2);
          if (gradientNorm <= m_options.gradientTolerance)
          {
            report.converged = true;
            report.reason = "gradient tolerance reached";
            break;
          }

          const auto oldNodes = m_geometry.nodes;
          Real step = m_options.initialStepSize;
          bool accepted = false;
          while (step >= m_options.minStepSize)
          {
            for (Index n = 0; n < m_geometry.nodes.size(); ++n)
            {
              if (m_geometry.nodes[n].fixed)
                continue;
              m_geometry.nodes[n].x[0] = oldNodes[n].x[0] - step * gradient[n][0];
              m_geometry.nodes[n].x[1] = oldNodes[n].x[1] - step * gradient[n][1];
            }

            const Real trial = m_objective.value(m_geometry);
            if (m_objective.invalidElementCount(m_geometry) == 0 &&
                std::isfinite(trial) &&
                trial < current)
            {
              current = trial;
              accepted = true;
              break;
            }

            m_geometry.nodes = oldNodes;
            step /= Real(2);
          }

          report.iterations = it + 1;
          if (!accepted)
          {
            report.reason = "line search failed";
            break;
          }
        }

        report.finalObjective = m_objective.value(m_geometry);
        report.minJacobian = m_objective.minJacobian(m_geometry);
        report.invalidElements = m_objective.invalidElementCount(m_geometry);
        if (report.reason.empty())
          report.reason = report.finalObjective < report.initialObjective
            ? "maximum iterations reached after objective reduction"
            : "maximum iterations reached";
        return report;
      }

    private:
      HighOrderTriangleGeometry& m_geometry;
      const Objective& m_objective;
      OptimizerOptions m_options;
  };
}

#endif
