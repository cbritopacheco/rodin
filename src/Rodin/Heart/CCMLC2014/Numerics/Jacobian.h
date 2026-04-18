/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jacobian.h
 * @brief Jacobian assembly for the CCMLC2014 coupled 0D residual.
 */
#ifndef RODIN_HEART_CCMLC2014_NUMERICS_JACOBIAN_H
#define RODIN_HEART_CCMLC2014_NUMERICS_JACOBIAN_H

#include <numbers>

#include "Rodin/Heart/CCMLC2014/Model/State.h"

namespace Rodin::Heart::CCMLC2014::Numerics
{
  /**
   * @brief Assemble the Jacobian of the fully coupled 0D nonlinear residual.
   */
  template <class Input, class DenseMatrix, class EvalData>
  inline void evaluateDynamicJacobian(
      const Input& in,
      DenseMatrix& jacobianMatrix,
      const EvalData& d,
      typename DenseMatrix::Scalar dt)
  {
    using Scalar = typename DenseMatrix::Scalar;
    jacobianMatrix.resize(Model::NVAR, Model::NVAR);
    jacobianMatrix.setZero();

    const Scalar inertiaCoefficient = in.d0 * in.rho;
    const Scalar geometricStretch = Scalar(1) + d.yMid / in.R0;
    const Scalar geometricStretchSquared = geometricStretch * geometricStretch;

    const Scalar totalStress =
      d.stressPassive + d.stressViscous + d.active.activeStress;
    const Scalar totalDiffStress =
      d.diffStressPassive + d.diffStressViscous + d.active.dActiveStressWrtDisplacement;

    jacobianMatrix(Model::DISP, Model::DISP) +=
      inertiaCoefficient / (dt * dt)
      + Scalar(0.5) * in.d0 / (in.R0 * in.R0) * totalStress
      + Scalar(0.5) * in.d0 / in.R0 * geometricStretch * totalDiffStress
      - Scalar(1) / in.R0 * d.pvMid * geometricStretch;
    jacobianMatrix(Model::DISP, Model::PV) += -Scalar(0.5) * geometricStretchSquared;

    jacobianMatrix(Model::PV, Model::DISP) +=
      Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * geometricStretch * d.vel
      + Scalar(4) * std::numbers::pi_v<Scalar> * in.R0 * in.R0 * geometricStretchSquared * (Scalar(1) / dt);

    jacobianMatrix(Model::PV, Model::PV) += in.cavityCapacity / dt + Scalar(0.5) * d.dCavityFluxCur_dPv;
    jacobianMatrix(Model::PV, Model::PAR) += Scalar(0.5) * d.dCavityFluxCur_dPar;

    jacobianMatrix(Model::PAR, Model::PV) += -d.dWindkesselOutflow_dPv;
    jacobianMatrix(Model::PAR, Model::PAR) +=
      in.Cp / dt + Scalar(1) / (Scalar(2) * in.Rp) - d.dWindkesselOutflow_dPar;
    jacobianMatrix(Model::PAR, Model::PD) += -Scalar(1) / (Scalar(2) * in.Rp);

    jacobianMatrix(Model::PD, Model::PAR) += -Scalar(1) / (Scalar(2) * in.Rp);
    jacobianMatrix(Model::PD, Model::PD) +=
      in.Cd / dt + Scalar(1) / (Scalar(2) * in.Rp) + Scalar(1) / (Scalar(2) * in.Rd);
  }
}

#endif
