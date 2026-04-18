/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Residual.h
 * @brief Residual assembly and intermediate data construction for CCMLC2014 0D dynamics.
 */
#ifndef RODIN_HEART_CCMLC2014_NUMERICS_RESIDUAL_H
#define RODIN_HEART_CCMLC2014_NUMERICS_RESIDUAL_H

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
   * @brief Build intermediate quantities used by residual and Jacobian assembly.
   */
  template <class PassiveLaw, class Input, class DenseVector, class State, class EvalData>
  inline void buildEvalData(
      const Input& in,
      const DenseVector& candidateUnknowns,
      const State& currentState,
      const State& previousState,
      typename DenseVector::Scalar tnp1,
      typename DenseVector::Scalar dt,
      EvalData& evalData)
  {
    using Scalar = typename DenseVector::Scalar;
    evalData.sn = currentState;
    evalData.snm1 = previousState;

    evalData.tnp1 = tnp1;
    evalData.dt = dt;

    evalData.y = candidateUnknowns[Model::DISP];
    evalData.pv = candidateUnknowns[Model::PV];
    evalData.par = candidateUnknowns[Model::PAR];
    evalData.pd = candidateUnknowns[Model::PD];

    evalData.yPrev = currentState.y;
    evalData.pvPrev = currentState.pv;
    evalData.parPrev = currentState.par;
    evalData.pdPrev = currentState.pd;
    evalData.yPrevPrev = previousState.y;

    evalData.pAtCur = in.pAt(tnp1);
    evalData.pAtPrev = in.pAt(currentState.t);
    evalData.pSvMid = in.pSv(currentState.t + Scalar(0.5) * dt);

    Physics::computeMidpointKinematics(in, evalData);
    Physics::computePassiveContribution<PassiveLaw>(in, evalData);
    Physics::computeViscousContribution(in, evalData);
    Physics::solveLocalDynamicActive(in, evalData, evalData.active);
    Physics::computeCavityFluxes(in, evalData);
    Physics::computeWindkesselOutflow(in, evalData);
  }

  /**
   * @brief Evaluate the coupled 0D residual vector at the current nonlinear iterate.
   */
  template <class Input, class DenseVector, class EvalData>
  inline void evaluateDynamicResidual(
      const Input& in,
      const EvalData& evalData,
      DenseVector& residualVector)
  {
    using Scalar = typename DenseVector::Scalar;
    residualVector.resize(Model::NVAR);
    residualVector.setZero();

    const Scalar inertiaCoefficient = in.d0 * in.rho;
    const Scalar geometricStretch = Scalar(1) + evalData.yMid / in.R0;
    const Scalar geometricStretchSquared = geometricStretch * geometricStretch;

    const Scalar totalStress =
      evalData.stressPassive + evalData.stressViscous + evalData.active.activeStress;

    residualVector[Model::DISP] =
      inertiaCoefficient / (evalData.dt * evalData.dt) * (evalData.y - Scalar(2) * evalData.yPrev + evalData.yPrevPrev)
      + in.d0 / in.R0 * geometricStretch * totalStress
      - evalData.pvMid * geometricStretchSquared;

    {
      const Scalar cavityCapacityCurrent = in.cavityCapacity / evalData.dt * evalData.pv;
      const Scalar cavityCapacityPrevious = in.cavityCapacity / evalData.dt * evalData.pvPrev;
      const Scalar volumeTerm =
        Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * in.R0 * geometricStretchSquared * evalData.vel;

      residualVector[Model::PV] =
        (cavityCapacityCurrent - cavityCapacityPrevious)
        + Scalar(0.5) * (evalData.cavityFluxCur + evalData.cavityFluxPrev)
        + volumeTerm;
    }

    {
      residualVector[Model::PAR] =
        in.Cp / evalData.dt * (evalData.par - evalData.parPrev)
        + (evalData.parMid - evalData.pdMid) / in.Rp
        - evalData.windkesselOutflow;
    }

    {
      residualVector[Model::PD] =
        in.Cd / evalData.dt * (evalData.pd - evalData.pdPrev)
        + (evalData.pdMid - evalData.parMid) / in.Rp
        - (evalData.pSvMid - evalData.pdMid) / in.Rd;
    }
  }
}

#endif
