/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DynamicSystem.h
 * @brief Residual and Jacobian assembly for the CCMLC2014 coupled 0D system.
 *
 * Provides the DynamicSystem class which orchestrates all physics evaluators
 * (kinematics, passive stress, viscous stress, active dynamics, cavity flux,
 * Windkessel outflow) and assembles the coupled nonlinear residual and
 * Jacobian matrices for the 0D left-ventricle model described in
 * Caruel et al. (2014), §4.
 *
 * The four coupled unknowns are radial displacement @f$ y @f$,
 * ventricular pressure @f$ p_v @f$, proximal arterial pressure @f$ p_{ar} @f$,
 * and distal pressure @f$ p_d @f$, as defined by Model::Variable.
 */
#ifndef RODIN_HEART_CCMLC2014_NUMERICS_DYNAMICSYSTEM_H
#define RODIN_HEART_CCMLC2014_NUMERICS_DYNAMICSYSTEM_H

#include <numbers>

#include "Rodin/Heart/CCMLC2014/Model/State.h"
#include "Rodin/Heart/CCMLC2014/Physics/ActiveDynamics.h"
#include "Rodin/Heart/CCMLC2014/Physics/Kinematics.h"
#include "Rodin/Heart/CCMLC2014/Physics/PassiveStress.h"
#include "Rodin/Heart/CCMLC2014/Physics/ValveLaw.h"
#include "Rodin/Heart/CCMLC2014/Physics/Windkessel.h"

namespace Rodin::Heart::CCMLC2014::Numerics
{
  /**
   * @brief Assembles the coupled residual and Jacobian for the 0D CCMLC2014 model.
   *
   * This class composes the individual physics evaluators and provides three
   * main operations:
   * - **buildEvalData**: populate the intermediate evaluation data structure
   *   from the candidate unknowns and time-history states.
   * - **evaluateResidual**: assemble the nonlinear residual vector.
   * - **evaluateJacobian**: assemble the consistent tangent (Jacobian) matrix.
   *
   * @tparam PassiveLaw Passive stress operator type (functor).
   * @tparam Input Model input parameter type.
   */
  template <class PassiveLaw, class Input>
  class DynamicSystem
  {
    public:
      /**
       * @brief Construct the dynamic system assembler.
       * @param[in] input Reference to model input parameters. Must outlive
       *   this object.
       */
      explicit DynamicSystem(const Input& input)
        : m_input(input),
          m_kinematics(input),
          m_passiveStress(input),
          m_viscousStress(input),
          m_activeDynamics(input),
          m_cavityFlux(input),
          m_windkessel(input)
      {}

      /**
       * @brief Build intermediate evaluation data from candidate unknowns.
       *
       * Populates @p evalData with midpoint kinematics, stress contributions,
       * active dynamics solution, cavity fluxes and Windkessel outflow.
       *
       * @tparam DenseVector Dense vector type for the unknown vector.
       * @tparam State State type for time-level data.
       * @tparam EvalData Evaluation data structure type.
       * @param[in] candidateUnknowns Current nonlinear iterate
       *   @f$ (y, p_v, p_{ar}, p_d)^{n+1} @f$.
       * @param[in] currentState State at time @f$ t_n @f$.
       * @param[in] previousState State at time @f$ t_{n-1} @f$.
       * @param[in] tnp1 Target time @f$ t_{n+1} @f$.
       * @param[in] dt Time-step size @f$ \Delta t @f$.
       * @param[out] evalData Populated evaluation data.
       */
      template <class DenseVector, class StateType, class EvalData>
      void buildEvalData(
          const DenseVector& candidateUnknowns,
          const StateType& currentState,
          const StateType& previousState,
          typename DenseVector::Scalar tnp1,
          typename DenseVector::Scalar dt,
          EvalData& evalData) const
      {
        using Scalar = typename DenseVector::Scalar;
        evalData.sn = currentState;
        evalData.snm1 = previousState;

        evalData.tnp1 = tnp1;
        evalData.dt = dt;

        evalData.y = candidateUnknowns[Model::RadialDisplacement];
        evalData.pv = candidateUnknowns[Model::VentricularPressure];
        evalData.par = candidateUnknowns[Model::ArterialPressure];
        evalData.pd = candidateUnknowns[Model::DistalPressure];

        evalData.yPrev = currentState.y;
        evalData.pvPrev = currentState.pv;
        evalData.parPrev = currentState.par;
        evalData.pdPrev = currentState.pd;
        evalData.yPrevPrev = previousState.y;

        evalData.pAtCur = m_input.pAt(tnp1);
        evalData.pAtPrev = m_input.pAt(currentState.t);
        evalData.pSvMid = m_input.pSv(currentState.t + Scalar(0.5) * dt);

        m_kinematics.evaluate(evalData);
        m_passiveStress.evaluate(evalData);
        m_viscousStress.evaluate(evalData);
        m_activeDynamics.evaluate(evalData, evalData.active);
        m_cavityFlux.evaluate(evalData);
        m_windkessel.evaluate(evalData);
      }

      /**
       * @brief Evaluate the coupled 0D residual vector.
       *
       * The residual encodes:
       * - **Row 0 (RadialDisplacement)**: wall momentum balance
       *   (inertia + stress divergence - pressure loading).
       * - **Row 1 (VentricularPressure)**: cavity mass balance
       *   (compliance + valve flux + volume rate).
       * - **Row 2 (ArterialPressure)**: proximal Windkessel balance
       *   (compliance + resistive flow - outflow).
       * - **Row 3 (DistalPressure)**: distal Windkessel balance.
       *
       * @tparam DenseVector Dense vector type.
       * @tparam EvalData Evaluation data structure type.
       * @param[in] evalData Fully populated evaluation data.
       * @param[out] residualVector Residual vector (resized to NumberOfVariables).
       */
      template <class DenseVector, class EvalData>
      void evaluateResidual(
          const EvalData& evalData,
          DenseVector& residualVector) const
      {
        using Scalar = typename DenseVector::Scalar;
        residualVector.resize(Model::NumberOfVariables);
        residualVector.setZero();

        const Scalar inertiaCoefficient = m_input.d0 * m_input.rho;
        const Scalar geometricStretch =
          Scalar(1) + evalData.yMid / m_input.R0;
        const Scalar geometricStretchSquared =
          geometricStretch * geometricStretch;

        const Scalar totalStress =
          evalData.stressPassive + evalData.stressViscous
          + evalData.active.activeStress;

        // Wall momentum balance
        residualVector[Model::RadialDisplacement] =
          inertiaCoefficient / (evalData.dt * evalData.dt)
            * (evalData.y - Scalar(2) * evalData.yPrev + evalData.yPrevPrev)
          + m_input.d0 / m_input.R0 * geometricStretch * totalStress
          - evalData.pvMid * geometricStretchSquared;

        // Cavity mass balance
        {
          const Scalar cavityCapacityCurrent =
            m_input.cavityCapacity / evalData.dt * evalData.pv;
          const Scalar cavityCapacityPrevious =
            m_input.cavityCapacity / evalData.dt * evalData.pvPrev;
          const Scalar volumeTerm =
            Scalar(4) * std::numbers::pi_v<Scalar>
            * m_input.R0 * m_input.R0
            * geometricStretchSquared * evalData.vel;

          residualVector[Model::VentricularPressure] =
            (cavityCapacityCurrent - cavityCapacityPrevious)
            + Scalar(0.5) * (evalData.cavityFluxCur + evalData.cavityFluxPrev)
            + volumeTerm;
        }

        // Proximal Windkessel balance
        residualVector[Model::ArterialPressure] =
          m_input.Cp / evalData.dt * (evalData.par - evalData.parPrev)
          + evalData.windkesselflowP;

        // Distal Windkessel balance
        residualVector[Model::DistalPressure] =
          m_input.Cd / evalData.dt * (evalData.pd - evalData.pdPrev)
          + evalData.windkesselflowD;
      }

      /**
       * @brief Assemble the Jacobian of the coupled 0D nonlinear residual.
       *
       * The Jacobian matrix @f$ J_{ij} = \partial R_i / \partial x_j @f$
       * is computed analytically from the intermediate evaluation data.
       *
       * @tparam DenseMatrix Dense matrix type.
       * @tparam EvalData Evaluation data structure type.
       * @param[in] evalData Fully populated evaluation data.
       * @param[out] jacobianMatrix Jacobian matrix (resized to
       *   NumberOfVariables × NumberOfVariables).
       * @param[in] dt Time-step size.
       */
      template <class DenseMatrix, class EvalData>
      void evaluateJacobian(
          const EvalData& evalData,
          DenseMatrix& jacobianMatrix,
          typename DenseMatrix::Scalar dt) const
      {
        using Scalar = typename DenseMatrix::Scalar;
        jacobianMatrix.resize(Model::NumberOfVariables, Model::NumberOfVariables);
        jacobianMatrix.setZero();

        const Scalar inertiaCoefficient = m_input.d0 * m_input.rho;
        const Scalar geometricStretch =
          Scalar(1) + evalData.yMid / m_input.R0;
        const Scalar geometricStretchSquared =
          geometricStretch * geometricStretch;

        const Scalar totalStress =
          evalData.stressPassive + evalData.stressViscous
          + evalData.active.activeStress;
        const Scalar totalDiffStress =
          evalData.diffStressPassive + evalData.diffStressViscous
          + evalData.active.dActiveStressWrtDisplacement;

        // --- Row: RadialDisplacement ---
        jacobianMatrix(Model::RadialDisplacement, Model::RadialDisplacement) +=
          inertiaCoefficient / (dt * dt)
          + Scalar(0.5) * m_input.d0 / (m_input.R0 * m_input.R0) * totalStress
          + Scalar(0.5) * m_input.d0 / m_input.R0
            * geometricStretch * totalDiffStress
          - Scalar(1) / m_input.R0 * evalData.pvMid * geometricStretch;

        jacobianMatrix(Model::RadialDisplacement, Model::VentricularPressure) +=
          -Scalar(0.5) * geometricStretchSquared;

        // --- Row: VentricularPressure ---
        jacobianMatrix(Model::VentricularPressure, Model::RadialDisplacement) +=
          Scalar(4) * std::numbers::pi_v<Scalar> * m_input.R0
            * geometricStretch * evalData.vel
          + Scalar(4) * std::numbers::pi_v<Scalar>
            * m_input.R0 * m_input.R0
            * geometricStretchSquared * (Scalar(1) / dt);

        jacobianMatrix(Model::VentricularPressure, Model::VentricularPressure) +=
          m_input.cavityCapacity / dt
          + Scalar(0.5) * evalData.dCavityFluxCur_dPv;

        jacobianMatrix(Model::VentricularPressure, Model::ArterialPressure) +=
          Scalar(0.5) * evalData.dCavityFluxCur_dPar;

        // --- Row: ArterialPressure ---
        jacobianMatrix(Model::ArterialPressure, Model::VentricularPressure) +=
          -evalData.dWindkesselOutflow_dPv;

        jacobianMatrix(Model::ArterialPressure, Model::ArterialPressure) +=
          m_input.Cp / dt + evalData.dWindkesselflowP_dPar;

        jacobianMatrix(Model::ArterialPressure, Model::DistalPressure) +=
          evalData.dWindkesselflowP_dPd;

        // --- Row: DistalPressure ---
        jacobianMatrix(Model::DistalPressure, Model::ArterialPressure) +=
          evalData.dWindkesselflowD_dPar;

        jacobianMatrix(Model::DistalPressure, Model::DistalPressure) +=
          m_input.Cd / dt
          + evalData.dWindkesselflowD_dPd;
      }

    private:
      const Input& m_input;
      Physics::MidpointKinematics<Input> m_kinematics;
      Physics::PassiveStressEvaluator<PassiveLaw, Input> m_passiveStress;
      Physics::ViscousStressEvaluator<Input> m_viscousStress;
      Physics::ActiveDynamicsSolver<Input> m_activeDynamics;
      Physics::CavityFluxEvaluator<Input> m_cavityFlux;
      Physics::WindkesselOutflowEvaluator<Input> m_windkessel;
  };
}

#endif
