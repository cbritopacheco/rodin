/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_QUALITYSAMPLES_H
#define RODIN_ADAPTATION_SWIFT_QUALITYSAMPLES_H

#include <algorithm>
#include <vector>

#include "Rodin/Variational/IntegrationPoint.h"
#include "QualityLattice.h"
#include "Parameters.h"
#include "../CellDeformation.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Common witnesses for affine hinges and actual quality checks.
   *
   * Uniform reference lattices include vertices and boundary points. Hinges
   * and actual quality checks use the same points. Equal positive reference
   * weights determine mapped cell mass; adaptive hinge weights redistribute
   * that mass. These witnesses are not degree-exact integration rules.
   */
  class QualitySamples
  {
    public:
      /**
       * @brief Constructs the cell quality sampling policy.
       * @param cell Reference cell.
       * @param order Displacement finite-element order.
       * @param parameters Quality sampling policy.
       */
      QualitySamples(const Geometry::Polytope& cell, size_t order,
        const Parameters& parameters)
        : m_cell(cell),
          m_parameters(parameters),
          m_formula(QualityLattice::get(cell.getGeometry(),
            parameters.quadrature.getQualityOrder(order,
              cell.getTransformation().getOrder(),
              Geometry::Polytope::Traits(cell.getGeometry()).getVertexCount() ==
                cell.getDimension() + 1)))
      {}

      /**
       * @brief Visits each quality witness and its positive penalty weight.
       * @param evaluate Callable accepting an IntegrationPoint and a weight.
       */
      template <class Evaluate>
      void forEach(Evaluate&& evaluate) const
      {
        const auto& cell = m_cell.get();
        const auto& formula = m_formula.get();
        const auto& quadrature = cell.getQuadrature(formula);
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          evaluate(Variational::IntegrationPoint(point, &formula, q),
            formula.getWeight(q) * point.getDistortion());
        }
      }

      /**
       * @brief Visits the frozen, constraint-specific adaptive hinge measures.
       * @param currentJacobian Frozen outer displacement gradient.
       * @param predictorJacobian Frozen directionally scaled predictor gradient.
       * @param evaluate Callable accepting a point, Jacobian weight and distortion weight.
       *
       * Each measure mixes equal weights with normalized guard penetration at
       * the current and full-predictor geometries. Both are strictly positive
       * and preserve mapped cell mass. No weight depends on the inner iterate.
       */
      template <class Jacobian, class Evaluate>
      void forEachHinge(const Jacobian& currentJacobian,
        const Jacobian& predictorJacobian, Evaluate&& evaluate) const
      {
        const auto& model = m_parameters.get().model;
        const size_t count = m_formula.get().getSize();
        std::vector<Real> riskJ(count), riskQ(count);
        Real measure = 0, sumJ = 0, sumQ = 0;
        size_t index = 0;
        CellDeformation state(m_cell.get().getDimension());
        const Real deltaJ = model.qualityGuard * (Real(1) - model.jacobian);
        const Real deltaQ = model.qualityGuard * (model.distortion - Real(1));
        forEach([&](const Variational::IntegrationPoint& ip, Real weight) {
          const Math::SpatialMatrix<Real> current = currentJacobian.getValue(ip);
          const Math::SpatialMatrix<Real> predictor = predictorJacobian.getValue(ip);
          Real jRisk = 0, qRisk = 0;
          for (const Real alpha : {Real(0), Real(1)})
          {
            state.setDisplacementGradient(
              Math::SpatialMatrix<Real>(current + alpha * predictor));
            jRisk = std::max(jRisk, std::clamp(
              Real(1) - (state.getJacobian() - model.jacobian) / deltaJ,
              Real(0), MaximumRisk));
            // Distortion is undefined at inverted predictor witnesses.
            qRisk = std::max(qRisk, state.isAdmissible()
              ? std::clamp(Real(1) -
                  (model.distortion - state.getRelativeDistortion()) / deltaQ,
                  Real(0), MaximumRisk)
              : MaximumRisk);
          }
          riskJ[index] = jRisk;
          riskQ[index] = qRisk;
          sumJ += jRisk;
          sumQ += qRisk;
          measure += weight;
          ++index;
        });
        const Real equalWeight = measure / count;
        const auto weight = [&](Real risk, Real sum) {
          return sum > Real(0)
            ? equalWeight * (Real(1) - AdaptiveFraction +
                AdaptiveFraction * count * risk / sum)
            : equalWeight;
        };
        index = 0;
        forEach([&](const Variational::IntegrationPoint& ip, Real) {
          evaluate(ip, weight(riskJ[index], sumJ), weight(riskQ[index], sumQ));
          ++index;
        });
      }

    private:
      /// Equal-weight fraction protects witnesses with no predicted risk.
      static constexpr Real AdaptiveFraction = Real(0.5);
      /// Caps dimensionless risk, including inverted predictor witnesses.
      static constexpr Real MaximumRisk = Real(100);
      std::reference_wrapper<const Geometry::Polytope> m_cell;
      std::reference_wrapper<const Parameters> m_parameters;
      std::reference_wrapper<const QualityLattice> m_formula;
  };
}

#endif
