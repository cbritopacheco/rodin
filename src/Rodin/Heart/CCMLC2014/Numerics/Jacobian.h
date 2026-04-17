/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_NUMERICS_JACOBIAN_H
#define RODIN_HEART_CCMLC2014_NUMERICS_JACOBIAN_H

#include <numbers>

#include "Rodin/Heart/CCMLC2014/Model/State.h"

namespace Rodin::Heart::CCMLC2014::Numerics
{
  template <class Input, class DenseMatrix, class EvalData>
  inline void evaluateDynamicJacobian(
      const Input& in,
      DenseMatrix& J,
      const EvalData& d,
      typename DenseMatrix::Scalar dt)
  {
    using Scalar = typename DenseMatrix::Scalar;
    J.resize(Model::NVAR, Model::NVAR);
    J.setZero();

    const Scalar coeff = in.d0 * in.rho;
    const Scalar geom = Scalar(1) + d.yMid / in.R0;
    const Scalar geom2 = geom * geom;

    const Scalar totalStress =
      d.stressPassive + d.stressViscous + d.active.stressActive;
    const Scalar totalDiffStress =
      d.diffStressPassive + d.diffStressViscous + d.active.diffStressActive;

    J(Model::DISP, Model::DISP) +=
      coeff / (dt * dt)
      + Scalar(0.5) * in.d0 / (in.R0 * in.R0) * totalStress
      + Scalar(0.5) * in.d0 / in.R0 * geom * totalDiffStress
      - Scalar(1) / in.R0 * d.pvMid * geom;
    J(Model::DISP, Model::PV) += -Scalar(0.5) * geom2;

    J(Model::PV, Model::DISP) +=
      Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * geom * d.vel
      + Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * in.R0 * geom2 * (Scalar(1) / dt);

    J(Model::PV, Model::PV) += in.cavityCapacity / dt + Scalar(0.5) * d.dCavityFluxCur_dPv;
    J(Model::PV, Model::PAR) += Scalar(0.5) * d.dCavityFluxCur_dPar;

    J(Model::PAR, Model::PV) += -d.dWindkesselOutflow_dPv;
    J(Model::PAR, Model::PAR) += in.Cp / dt + Scalar(1) / (Scalar(2) * in.Rp) - d.dWindkesselOutflow_dPar;
    J(Model::PAR, Model::PD) += -Scalar(1) / (Scalar(2) * in.Rp);

    J(Model::PD, Model::PAR) += -Scalar(1) / (Scalar(2) * in.Rp);
    J(Model::PD, Model::PD) += in.Cd / dt + Scalar(1) / (Scalar(2) * in.Rp) + Scalar(1) / (Scalar(2) * in.Rd);
  }
}

#endif
