/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_OBJECTIVE_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_OBJECTIVE_H

#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"

#include "Geometry.h"
#include "Metrics.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  class InterfaceFitPenalty
  {
    public:
      virtual ~InterfaceFitPenalty() = default;
      virtual Real value(const HighOrderTriangleGeometry&) const
      {
        return Real(0);
      }
  };

  class Objective
  {
    public:
      /**
       * @brief Builds a fixed-topology objective from a quality metric.
       *
       * The current objective is:
       *   quality(geometry) + mu * deviation_from_initial_geometry.
       *
       * Rodin-native TMOP terms expose the corresponding residual and tangent
       * contributions for LinearProblem assembly.
       */
      explicit Objective(const MetricBase& metric)
        : m_metric(metric)
      {}

      /**
       * @brief Sets the weight of the quadratic deviation penalty.
       */
      Objective& setDeviationWeight(Real weight)
      {
        m_deviationWeight = std::max(Real(0), weight);
        return *this;
      }

      /**
       * @brief Sets the reference sample points used for objective evaluation.
       */
      Objective& setSamplePoints(std::vector<ReferencePoint> points)
      {
        if (points.empty())
          Alert::MemberFunctionException(*this, __func__)
            << "Expected at least one TMOP sample point."
            << Alert::Raise;
        m_points = std::move(points);
        return *this;
      }

      /**
       * @brief Evaluates quality plus deviation.
       */
      Real value(const HighOrderTriangleGeometry& geometry) const
      {
        return qualityValue(geometry) + deviationValue(geometry);
      }

      /**
       * @brief Evaluates the quality term only.
       */
      Real qualityValue(const HighOrderTriangleGeometry& geometry) const
      {
        Real res = 0;
        CurvedTriangleJacobianEvaluator evaluator;
        for (Index c = 0; c < geometry.cells.size(); ++c)
          for (const auto& point : m_points)
            res += m_metric.value(evaluator.jacobian(geometry, c, point));
        return res / static_cast<Real>(m_points.size());
      }

      /**
       * @brief Evaluates the quadratic deviation penalty only.
       */
      Real deviationValue(const HighOrderTriangleGeometry& geometry) const
      {
        if (m_deviationWeight == Real(0))
          return Real(0);

        Real res = 0;
        for (const auto& node : geometry.nodes)
        {
          const auto dx = node.x - node.x0;
          res += dx.dot(dx);
        }
        return m_deviationWeight * res;
      }

      /**
       * @brief Minimum sampled geometry Jacobian determinant.
       */
      Real minJacobian(const HighOrderTriangleGeometry& geometry) const
      {
        return CurvedTriangleJacobianEvaluator().minJacobian(geometry, m_points);
      }

      /**
       * @brief Number of cells with at least one invalid sampled Jacobian.
       */
      Index invalidElementCount(const HighOrderTriangleGeometry& geometry) const
      {
        Index res = 0;
        CurvedTriangleJacobianEvaluator evaluator;
        for (Index c = 0; c < geometry.cells.size(); ++c)
        {
          bool invalid = false;
          for (const auto& point : m_points)
            invalid = invalid ||
              evaluator.determinant(geometry, c, point) <=
                m_metric.getDeterminantTolerance();
          if (invalid)
            res++;
        }
        return res;
      }

    private:
      const MetricBase& m_metric;
      Real m_deviationWeight = 0;
      std::vector<ReferencePoint> m_points = {{
        { Real(1) / Real(3), Real(1) / Real(3) },
        { Real(0.2), Real(0.2) },
        { Real(0.6), Real(0.2) },
        { Real(0.2), Real(0.6) } }};
  };
}

#endif
