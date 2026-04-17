/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
  template <class PassiveLaw, class Input, class DenseVector, class State, class EvalData>
  inline void buildEvalData(
      const Input& in,
      const DenseVector& x,
      const State& sn,
      const State& snm1,
      typename DenseVector::Scalar tnp1,
      typename DenseVector::Scalar dt,
      EvalData& d)
  {
    using Scalar = typename DenseVector::Scalar;
    d.sn = sn;
    d.snm1 = snm1;

    d.tnp1 = tnp1;
    d.dt = dt;

    d.y = x[Model::DISP];
    d.pv = x[Model::PV];
    d.par = x[Model::PAR];
    d.pd = x[Model::PD];

    d.yPrev = sn.y;
    d.pvPrev = sn.pv;
    d.parPrev = sn.par;
    d.pdPrev = sn.pd;
    d.yPrevPrev = snm1.y;

    d.pAtCur = in.pAt(tnp1);
    d.pAtPrev = in.pAt(sn.t);
    d.pSvMid = in.pSv(sn.t + Scalar(0.5) * dt);

    Physics::computeMidpointKinematics(in, d);
    Physics::computePassiveContribution<PassiveLaw>(in, d);
    Physics::computeViscousContribution(in, d);
    Physics::solveLocalDynamicActive(in, d, d.active);
    Physics::computeCavityFluxes(in, d);
    Physics::computeWindkesselOutflow(in, d);
  }

  template <class Input, class DenseVector, class EvalData>
  inline void evaluateDynamicResidual(
      const Input& in,
      const EvalData& d,
      DenseVector& R)
  {
    using Scalar = typename DenseVector::Scalar;
    R.resize(Model::NVAR);
    R.setZero();

    const Scalar coeff = in.d0 * in.rho;
    const Scalar geom = Scalar(1) + d.yMid / in.R0;
    const Scalar geom2 = geom * geom;

    const Scalar totalStress =
      d.stressPassive + d.stressViscous + d.active.stressActive;

    R[Model::DISP] =
      coeff / (d.dt * d.dt) * (d.y - Scalar(2) * d.yPrev + d.yPrevPrev)
      + in.d0 / in.R0 * geom * totalStress
      - d.pvMid * geom2;

    {
      const Scalar capacityCur = in.cavityCapacity / d.dt * d.pv;
      const Scalar capacityPrev = in.cavityCapacity / d.dt * d.pvPrev;
      const Scalar volumeTerm =
        Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * in.R0 * geom2 * d.vel;

      R[Model::PV] =
        (capacityCur - capacityPrev)
        + Scalar(0.5) * (d.cavityFluxCur + d.cavityFluxPrev)
        + volumeTerm;
    }

    {
      R[Model::PAR] =
        in.Cp / d.dt * (d.par - d.parPrev)
        + (d.parMid - d.pdMid) / in.Rp
        - d.windkesselOutflow;
    }

    {
      R[Model::PD] =
        in.Cd / d.dt * (d.pd - d.pdPrev)
        + (d.pdMid - d.parMid) / in.Rp
        - (d.pSvMid - d.pdMid) / in.Rd;
    }
  }
}

#endif
