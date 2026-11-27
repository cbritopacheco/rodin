/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_HINGE_WEIGHTS_H
#define RODIN_ADAPTATION_SWIFT_HINGE_WEIGHTS_H

#include <algorithm>
#include <vector>
#include "../../CellDeformation.h"
#include "../Witnesses.h"

namespace Rodin::Adaptation::SWIFT::Hinge
{
  /// @brief Mapped hinge measures at one witness.
  struct Measure
  {
    Real jacobian; ///< Weight of the Jacobian hinge.
    Real distortion; ///< Weight of the distortion hinge.
  };

  /// @brief Frozen per-cell measures, each constraint preserving the mapped cell mass.
  class Weights
  {
    public:
      /**
       * @brief Redistributes mapped cell mass without changing its total for either constraint.
       * @param witnesses Ordered geometric covering of the cell.
       * @param currentJacobian Displacement gradient at the frozen outer state.
       * @param predictorJacobian Gradient of the scaled predictor increment.
       * @param model Quality bounds and guard width.
       */
      template <class Jacobian>
      Weights(const Witnesses& witnesses, const Jacobian& currentJacobian,
        const Jacobian& predictorJacobian, const Parameters::Model& model)
      {
        const size_t count = witnesses.getSize();
        std::vector<Measure> risks(count);
        Real measure = 0, sumJ = 0, sumQ = 0;
        CellDeformation state(
          witnesses.getPoint(0).getPoint().getPolytope().getDimension());
        const Real deltaJ = model.qualityGuard * (Real(1) - model.jacobian);
        const Real deltaQ = model.qualityGuard * (model.distortion - Real(1));
        for (size_t q = 0; q < count; ++q)
        {
          const auto ip = witnesses.getPoint(q);
          const Math::SpatialMatrix<Real> current = currentJacobian.getValue(ip);
          const Math::SpatialMatrix<Real> predictor = predictorJacobian.getValue(ip);
          Real jRisk = 0, qRisk = 0;
          for (const Real alpha : {Real(0), Real(1)})
          {
            state.setDisplacementGradient(
              Math::SpatialMatrix<Real>(current + alpha * predictor));
            jRisk = std::max(jRisk,
              std::clamp(Real(1) - (state.getJacobian() - model.jacobian) / deltaJ,
                Real(0), MaximumRisk));
            qRisk = std::max(qRisk,
              state.isAdmissible()
                ? std::clamp(
                    Real(1) - (model.distortion - state.getRelativeDistortion()) / deltaQ,
                    Real(0), MaximumRisk)
                : MaximumRisk);
          }
          risks[q] = {jRisk, qRisk};
          sumJ += jRisk;
          sumQ += qRisk;
          measure += witnesses.getWeight(q);
        }
        m_uniform = measure / count;
        if (sumJ == Real(0) && sumQ == Real(0))
          return;
        const auto weight = [&](Real risk, Real sum) {
          return sum > Real(0) ? m_uniform *
              (Real(1) - AdaptiveFraction + AdaptiveFraction * count * risk / sum)
                               : m_uniform;
        };
        for (auto& risk : risks)
          risk = {weight(risk.jacobian, sumJ), weight(risk.distortion, sumQ)};
        m_measures = std::move(risks);
      }

      /// @returns Whether both constraints use one uniform weight.
      bool isUniform() const
      {
        return m_measures.empty();
      }

      /**
       * @param q Witness index in covering order.
       * @returns Both measures at the given witness index.
       */
      Measure operator[](size_t q) const
      {
        return isUniform() ? Measure{m_uniform, m_uniform} : m_measures[q];
      }

    private:
      static constexpr Real AdaptiveFraction = Real(0.5);
      static constexpr Real MaximumRisk = Real(100);
      Real m_uniform;
      std::vector<Measure> m_measures;
  };
}

#endif
