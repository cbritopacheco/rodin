/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PassiveStress.h
 * @brief Passive and viscous stress evaluators for the CCMLC2014 reduced 0D dynamics.
 *
 * Provides the PassiveStressEvaluator and ViscousStressEvaluator classes.
 * The passive stress is derived from the reduced strain-energy density
 * @f$ W(J_1, J_2, J_4) @f$ via chain rule through the reduced invariants
 * (see Caruel et al. 2014, §2.2 and §4). The viscous stress models
 * viscoelastic damping in the thin-wall shell.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_PASSIVESTRESS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_PASSIVESTRESS_H

#include <cmath>

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Evaluates the passive stress contribution from the constitutive law.
   *
   * Computes the passive second Piola-Kirchhoff-type stress and its tangent
   * (derivative with respect to the midpoint displacement) by evaluating the
   * passive energy law at the reduced right Cauchy-Green scalar @f$ C @f$.
   *
   * @tparam PassiveLaw Functor type implementing the passive stress operator.
   * @tparam Input Model input parameter type.
   */
  template <class PassiveLaw, class Input>
  class PassiveStressEvaluator
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters (passive energy law, geometry).
       */
      explicit PassiveStressEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate passive stress and its displacement derivative.
       *
       * Reads @p data.C, @p data.sqrtC and writes @p data.stressPassive,
       * @p data.diffStressPassive.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data with kinematics already computed.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar dC_dyMid = Scalar(2) * data.sqrtC / m_input.R0;

        Scalar passiveStress = 0.0;
        Scalar passiveStressDerivativeWrtMidpointDisplacement = 0.0;
        PassiveLaw passiveLaw;
        passiveLaw(
            m_input.passiveEnergy,
            data.C,
            dC_dyMid,
            passiveStress,
            passiveStressDerivativeWrtMidpointDisplacement);

        data.stressPassive = passiveStress;
        data.diffStressPassive =
          Scalar(0.5) * passiveStressDerivativeWrtMidpointDisplacement;
      }

    private:
      const Input& m_input;
  };

  /**
   * @brief Evaluates the viscous stress contribution.
   *
   * Computes viscous (rate-dependent) stress and its tangent with respect
   * to the midpoint displacement. The viscous model uses the time-rate of
   * the reduced Green-Lagrange strain components.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class ViscousStressEvaluator
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters (viscosity coefficient, geometry).
       */
      explicit ViscousStressEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate viscous stress and its displacement derivative.
       *
       * Reads kinematic fields from @p data and writes @p data.stressViscous,
       * @p data.diffStressViscous.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data with kinematics already computed.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar referenceRadius = m_input.R0;
        const Scalar viscosity = m_input.eta;
        const Scalar radialVelocity = data.vel;
        const Scalar reducedStretch = data.sqrtC;
        const Scalar reducedCauchyGreen = data.C;

        const Scalar dGreenRRWrtDisplacement =
          Scalar(-2) / referenceRadius * std::pow(reducedStretch, -5);
        const Scalar dGreenThetaThetaWrtDisplacement = data.diffGreen;

        const Scalar viscousStressRR =
          viscosity * dGreenRRWrtDisplacement * radialVelocity;
        const Scalar viscousStressThetaTheta =
          viscosity * dGreenThetaThetaWrtDisplacement * radialVelocity;
        const Scalar viscousStressPhiPhi = viscousStressThetaTheta;

        data.stressViscous =
          viscousStressThetaTheta + viscousStressPhiPhi
          - Scalar(2) * std::pow(reducedCauchyGreen, -3) * viscousStressRR;

        const Scalar dDotCRRWrtDisplacement =
          Scalar(10) / (referenceRadius * referenceRadius)
            * viscosity * std::pow(reducedStretch, -6) * radialVelocity
          + Scalar(2) * viscosity * dGreenRRWrtDisplacement / data.dt;

        const Scalar dDotCThetaThetaWrtDisplacement =
          viscosity / (referenceRadius * referenceRadius) * radialVelocity
          + Scalar(2) * viscosity * dGreenThetaThetaWrtDisplacement / data.dt;

        const Scalar viscousStressDerivativeWrtDisplacement =
          dDotCThetaThetaWrtDisplacement + dDotCThetaThetaWrtDisplacement
          - Scalar(2)
            * (std::pow(reducedCauchyGreen, -3) * dDotCRRWrtDisplacement
               - Scalar(6) / referenceRadius
                 * std::pow(reducedStretch, -7) * viscousStressRR);

        data.diffStressViscous = viscousStressDerivativeWrtDisplacement;
      }

    private:
      const Input& m_input;
  };
}

#endif
