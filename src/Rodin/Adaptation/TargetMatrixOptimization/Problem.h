/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_PROBLEM_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_PROBLEM_H

#include <functional>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Types.h"

#include "Geometry.h"
#include "Metrics.h"
#include "Objective.h"
#include "Optimizer.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  struct ProblemReport
  {
    OptimizationReport optimization;
    Real initialObjective = 0;
    Real finalObjective = 0;
    Real initialMinJacobian = 0;
    Real finalMinJacobian = 0;
    Index initialInvalidElements = 0;
    Index finalInvalidElements = 0;
  };

  /**
   * @brief Legacy standalone diagnostic facade for fixed-topology adaptation.
   *
   * This class predates the Rodin-native residual/tangent TMOP path. It reuses
   * Rodin mesh/connectivity infrastructure, constructs a detached order-2
   * geometry representation, configures the prototype objective, and advances
   * that detached geometry while keeping topology fixed. It is kept only as a
   * diagnostic bridge while the production path moves to QualityTerm,
   * DeviationTerm, LevelSetFitTerm, Variational::Problem, and NewtonSolver.
   *
   * Example:
   * @code
   * TMOP::SquaredDistanceMetric metric;
   * TMOP::Problem problem(mesh);
   * auto report = problem
   *   .setMetric(metric)
   *   .setGeometryOrder(2)
   *   .setDeviationWeight(1e-2)
   *   .initialize()
   *   .optimize();
   * @endcode
   *
   * The production API is not this class. Use the Rodin-native TMOP terms for
   * residual/tangent assembly into LinearProblem/Newton infrastructure.
   */
  class Problem
  {
    public:
      explicit Problem(const Geometry::LocalMesh& mesh)
        : m_mesh(mesh)
      {}

      Problem& setMetric(const MetricBase& metric)
      {
        m_metric = std::cref(metric);
        return *this;
      }

      Problem& setGeometryOrder(Index order)
      {
        m_geometryOrder = order;
        return *this;
      }

      Problem& setFixOriginalVertices(bool fixed)
      {
        m_fixOriginalVertices = fixed;
        return *this;
      }

      Problem& setDeviationWeight(Real weight)
      {
        m_deviationWeight = std::max(Real(0), weight);
        return *this;
      }

      Problem& setSamplePoints(std::vector<ReferencePoint> points)
      {
        if (points.empty())
          Alert::MemberFunctionException(*this, __func__)
            << "Expected at least one TMOP sample point."
            << Alert::Raise;
        m_samplePoints = std::move(points);
        return *this;
      }

      Problem& setOptimizerOptions(const OptimizerOptions& options)
      {
        m_optimizerOptions = options;
        return *this;
      }

      Problem& initialize()
      {
        m_geometry = HighOrderGeometryUpgrade()
          .setFixOriginalVertices(m_fixOriginalVertices)
          .upgrade(m_mesh.get(), m_geometryOrder);
        m_initialized = true;
        return *this;
      }

      bool isInitialized() const
      {
        return m_initialized;
      }

      HighOrderTriangleGeometry& getGeometry()
      {
        requireInitialized(__func__);
        return m_geometry;
      }

      const HighOrderTriangleGeometry& getGeometry() const
      {
        requireInitialized(__func__);
        return m_geometry;
      }

      Objective makeObjective() const
      {
        Objective objective(m_metric.get());
        objective
          .setDeviationWeight(m_deviationWeight)
          .setSamplePoints(m_samplePoints);
        return objective;
      }

      Real value() const
      {
        requireInitialized(__func__);
        return makeObjective().value(m_geometry);
      }

      Real minJacobian() const
      {
        requireInitialized(__func__);
        return makeObjective().minJacobian(m_geometry);
      }

      Index invalidElementCount() const
      {
        requireInitialized(__func__);
        return makeObjective().invalidElementCount(m_geometry);
      }

      ProblemReport optimize()
      {
        requireInitialized(__func__);

        ProblemReport report;
        auto objective = makeObjective();
        report.initialObjective = objective.value(m_geometry);
        report.initialMinJacobian = objective.minJacobian(m_geometry);
        report.initialInvalidElements = objective.invalidElementCount(m_geometry);

        report.optimization = Optimizer(m_geometry, objective)
          .setOptions(m_optimizerOptions)
          .optimize();

        report.finalObjective = objective.value(m_geometry);
        report.finalMinJacobian = objective.minJacobian(m_geometry);
        report.finalInvalidElements = objective.invalidElementCount(m_geometry);
        return report;
      }

    private:
      void requireInitialized(const char* function) const
      {
        if (!m_initialized)
          Alert::MemberFunctionException(*this, function)
            << "TMOP problem geometry has not been initialized."
            << Alert::Raise;
      }

      std::reference_wrapper<const Geometry::LocalMesh> m_mesh;
      SquaredDistanceMetric m_defaultMetric;
      std::reference_wrapper<const MetricBase> m_metric = m_defaultMetric;
      Index m_geometryOrder = 2;
      bool m_fixOriginalVertices = true;
      Real m_deviationWeight = 0;
      OptimizerOptions m_optimizerOptions;
      std::vector<ReferencePoint> m_samplePoints = {{
        { Real(1) / Real(3), Real(1) / Real(3) },
        { Real(0.2), Real(0.2) },
        { Real(0.6), Real(0.2) },
        { Real(0.2), Real(0.6) } }};
      HighOrderTriangleGeometry m_geometry;
      bool m_initialized = false;
  };
}

#endif
