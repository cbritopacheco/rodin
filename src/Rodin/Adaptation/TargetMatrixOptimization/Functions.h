/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_FUNCTIONS_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_FUNCTIONS_H

#include <functional>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/QuadratureRule.h"

#include "Metrics.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  /**
   * @brief Rodin FunctionBase adapter for TMOP metric evaluation.
   *
   * This object makes a TMOP metric usable as a normal Rodin integrand over a
   * mesh:
   *
   * @code
   * TMOP::SquaredDistanceMetric metric;
   * TMOP::JacobianMetricFunction w(metric);
   * Variational::QuadratureRule<Variational::FunctionBase<JacobianMetricFunction>> q(w);
   * q.setPolytope(cell);
   * auto local = q.compute();
   * @endcode
   *
   * It evaluates the metric on the Jacobian carried by Rodin's
   * Geometry::Point. This path is for linear/isoparametric mesh geometry and is
   * separate from the high-order geometry-node evaluator used by
   * TMOP::Objective. Both share the same metric value objects.
   */
  class JacobianMetricFunction final
    : public Variational::FunctionBase<JacobianMetricFunction>
  {
    public:
      explicit JacobianMetricFunction(const MetricBase& metric)
        : m_metric(metric)
      {}

      Real getValue(const Geometry::Point& point) const
      {
        const auto& J = point.getJacobian();
        Matrix2 A(2, 2);
        A(0, 0) = J(0, 0);
        A(0, 1) = J(0, 1);
        A(1, 0) = J(1, 0);
        A(1, 1) = J(1, 1);
        return m_metric.get().value(A);
      }

      Real getValue(const Variational::IntegrationPoint& ip) const
      {
        return getValue(ip.getPoint());
      }

      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      JacobianMetricFunction* copy() const noexcept override
      {
        return new JacobianMetricFunction(*this);
      }

    private:
      std::reference_wrapper<const MetricBase> m_metric;
  };

  /**
   * @brief Computes a TMOP metric objective through Rodin quadrature machinery.
   *
   * This is the minimal bridge to Rodin's FunctionBase/QuadratureRule
   * architecture. It is useful for linear mesh geometry, tests, and as a future
   * target for integrator-style nonlinear assembly.
   */
  class LinearMeshMetricObjective
  {
    public:
      explicit LinearMeshMetricObjective(const MetricBase& metric)
        : m_metric(metric)
      {}

      Real compute(const Geometry::LocalMesh& mesh) const
      {
        JacobianMetricFunction integrand(m_metric.get());
        Variational::QuadratureRule<
          Variational::FunctionBase<JacobianMetricFunction>> q(integrand);

        Real value = 0;
        for (auto it = mesh.getCell(); it; ++it)
        {
          const auto cell = *it;
          q.setPolytope(cell);
          value += q.compute();
        }
        return value;
      }

    private:
      std::reference_wrapper<const MetricBase> m_metric;
  };
}

#endif
