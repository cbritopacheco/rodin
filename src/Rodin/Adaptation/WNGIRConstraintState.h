/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRCONSTRAINTSTATE_H
#define RODIN_ADAPTATION_WNGIRCONSTRAINTSTATE_H

#include "CellDeformation.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Pointwise weights and affine targets of the WNGIR constraint model.
  class WNGIRConstraintState
  {
    public:
      /// @brief Constructs the w n g i r constraint state.
      WNGIRConstraintState(const CellDeformation& deformation,
        const WNGIRParameters& parameters,
        const Math::SpatialMatrix<Real>* predictorGradient = nullptr)
      {
        const Real j = deformation.getJacobian();
        const bool shapeOK = deformation.isAdmissible();
        Real predictedAJ = 0;
        Real predictedAQ = 0;
        if (predictorGradient && deformation.isInvertible())
        {
          predictedAJ = deformation.getJacobianAction(*predictorGradient);
          if (shapeOK)
            predictedAQ = deformation.getRelativeDistortionAction(*predictorGradient);
        }

        if (shapeOK)
        {
          const Real q = deformation.getRelativeDistortion();
          m_marginJ = j - parameters.jSafe;
          m_marginQ = parameters.qMax - q;
          if (parameters.includeAdmissibilityMetric && parameters.gammaJ > Real(0) &&
            m_marginJ > Real(0))
          {
            if (parameters.metricActivation == WNGIRMetricActivation::Smooth)
            {
              const Real s = std::max(m_marginJ, parameters.metricActivationEpsilon);
              m_jWeight =
                smoothStep(parameters.s0J - m_marginJ, parameters.jBarrierSmoothDelta) /
                (s * s);
            }
            else if (m_marginJ < parameters.s0J)
              m_jWeight = Real(1) / (m_marginJ * m_marginJ);
          }
          if (parameters.includeAdmissibilityMetric && parameters.gammaQ > Real(0) &&
            m_marginQ > Real(0))
          {
            if (parameters.metricActivation == WNGIRMetricActivation::Smooth)
            {
              const Real s = std::max(m_marginQ, parameters.metricActivationEpsilon);
              m_qWeight =
                smoothStep(parameters.s0Q - m_marginQ, parameters.qBarrierSmoothDelta) /
                (s * s);
            }
            else if (m_marginQ < parameters.s0Q)
              m_qWeight = Real(1) / (m_marginQ * m_marginQ);
          }
          if (parameters.includeQualityMetric && parameters.gammaQual > Real(0))
          {
            if (parameters.metricActivation == WNGIRMetricActivation::Smooth)
              m_qualityWeight =
                smoothStep(q - parameters.qStar, parameters.qualitySmoothDelta);
            else if (q > parameters.qStar)
              m_qualityWeight = Real(1);
          }

          const bool fractionalMargin = parameters.constraintFormulation ==
              WNGIRConstraintFormulation::FractionalMarginMetric ||
            parameters.constraintFormulation ==
              WNGIRConstraintFormulation::SafeguardedMarginMetric;
          const Real marginFraction =
            std::clamp(parameters.marginFraction, Real(0), Real(1));
          if (predictorGradient)
          {
            const Real jLimit = fractionalMargin ? marginFraction * m_marginJ : Real(0);
            const Real qLimit = fractionalMargin ? marginFraction * m_marginQ : Real(0);
            if (-predictedAJ <= jLimit)
              m_jWeight = 0;
            if (predictedAQ <= qLimit)
              m_qWeight = 0;
            if (predictedAQ <= Real(0))
              m_qualityWeight = 0;
          }
        }

        if (parameters.includeQualityMetric && parameters.gammaSize > Real(0))
        {
          if (parameters.metricActivation == WNGIRMetricActivation::Smooth)
            m_sizeWeight =
              smoothStep(parameters.jStar - j, parameters.jBarrierSmoothDelta);
          else if (j < parameters.jStar)
            m_sizeWeight = Real(1);
          if (predictorGradient && predictedAJ >= Real(0))
            m_sizeWeight = 0;
        }

        if (parameters.constraintFormulation ==
            WNGIRConstraintFormulation::FractionalMarginMetric ||
          parameters.constraintFormulation ==
            WNGIRConstraintFormulation::SafeguardedMarginMetric)
        {
          const Real fraction = std::clamp(parameters.marginFraction, Real(0), Real(1));
          m_jTarget = fraction * m_marginJ;
          m_qTarget = fraction * m_marginQ;
        }
      }

      /// @brief The jacobian weight.
      Real getJacobianWeight() const
      {
        return m_jWeight;
      }
      /// @brief The distortion weight.
      Real getDistortionWeight() const
      {
        return m_qWeight;
      }
      /// @brief The quality weight.
      Real getQualityWeight() const
      {
        return m_qualityWeight;
      }
      /// @brief The size weight.
      Real getSizeWeight() const
      {
        return m_sizeWeight;
      }
      /// @brief The jacobian target.
      Real getJacobianTarget() const
      {
        return m_jTarget;
      }
      /// @brief The distortion target.
      Real getDistortionTarget() const
      {
        return m_qTarget;
      }

    private:
      static Real smoothStep(Real x, Real delta)
      {
        if (delta <= Real(0))
          return x > Real(0) ? Real(1) : Real(0);
        return Real(0.5) * (Real(1) + std::tanh(x / delta));
      }

      Real m_marginJ = 0;
      Real m_marginQ = 0;
      Real m_jWeight = 0;
      Real m_qWeight = 0;
      Real m_qualityWeight = 0;
      Real m_sizeWeight = 0;
      Real m_jTarget = 0;
      Real m_qTarget = 0;
  };
}

#endif
