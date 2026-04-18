/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_MODEL_STATE_H
#define RODIN_HEART_CCMLC2014_MODEL_STATE_H

#include <cstddef>
#include <functional>

#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/PartialPivLU.h"

namespace Rodin::Heart::CCMLC2014::Model
{
  enum Variable : size_t
  {
    DISP = 0,
    PV,
    PAR,
    PD,
    NVAR
  };

  template <class Scalar>
  struct StateT
  {
    Scalar y = 0.0;
    Scalar v = 0.0;
    Scalar pv = 0.0;
    Scalar par = 0.0;
    Scalar pd = 0.0;

    Scalar ec = 0.0;
    Scalar kc = 0.0;
    Scalar tauc = 0.0;
    Scalar gamma = 0.0;
    Scalar beta = 0.0;

    Scalar t = 0.0;
  };

  template <class Scalar, class PassiveEnergyLaw>
  struct InputT
  {
    Scalar rho = 1.0;
    Scalar d0 = 1.0;
    Scalar R0 = 1.0;

    Scalar Es = 1.0;
    Scalar eta = 0.0;
    Scalar mu = 0.0;
    Scalar alpha = 0.0;
    Scalar k0 = 0.0;
    Scalar sigma0 = 0.0;

    Scalar Cp = 1.0;
    Scalar Cd = 1.0;
    Scalar Rp = 1.0;
    Scalar Rd = 1.0;

    Scalar Kat = 0.0;
    Scalar Kp = 0.0;
    Scalar Kar = 0.0;

    Scalar cavityCapacity = Scalar(5e-12);

    Scalar localTolerance = Scalar(1e-12);
    size_t localMaxIterations = 50;
    Scalar localDamping = Scalar(1.0);
    Scalar absRegularization = Scalar(1e-14);

    Scalar initFibDef = 0.0;
    Scalar initActiveStiffness = 0.0;
    Scalar initActiveStress = 0.0;

    std::function<Scalar(Scalar)> u =
      [](Scalar) { return Scalar(0); };

    std::function<Scalar(Scalar)> pAt =
      [](Scalar) { return Scalar(0); };

    std::function<Scalar(Scalar)> pSv =
      [](Scalar) { return Scalar(0); };

    PassiveEnergyLaw passiveEnergy;
  };

  template <class Scalar, class DenseLinearSystem>
  struct ReportT
  {
    bool converged = false;
    size_t iterations = 0;
    Scalar finalResidual = 0.0;
    Scalar finalStepNorm = 0.0;
    typename Solver::NewtonSolver<::Rodin::Solver::PartialPivLU<DenseLinearSystem>>::ConvergenceReason reason =
      Solver::NewtonSolver<::Rodin::Solver::PartialPivLU<DenseLinearSystem>>::ConvergenceReason::MaxIterations;
  };

  template <class Scalar>
  struct LocalActiveDataT
  {
    Scalar fib0 = 0.0;
    Scalar fib1 = 0.0;
    Scalar fib12 = 0.0;

    Scalar gammaOld = 0.0;
    Scalar betaOld = 0.0;
    Scalar gammaNew = 0.0;
    Scalar betaNew = 0.0;

    Scalar u1 = 0.0;
    Scalar u1Plus = 0.0;
    Scalar n0 = 0.0;

    Scalar k21 = 0.0;
    Scalar k22 = 0.0;
    Scalar krc_k22 = 0.0;

    Scalar sigma1d = 0.0;
    Scalar partialSigma1dWrtDisp = 0.0;
    Scalar partialSigma1dWrtEc = 0.0;

    Scalar stressActive = 0.0;
    Scalar diffStressActive = 0.0;

    bool converged = false;
    size_t iterations = 0;
  };

  template <class Scalar>
  struct EvalDataT
  {
    StateT<Scalar> sn;
    StateT<Scalar> snm1;

    Scalar tnp1 = 0.0;
    Scalar dt = 0.0;

    Scalar y = 0.0;
    Scalar pv = 0.0;
    Scalar par = 0.0;
    Scalar pd = 0.0;

    Scalar yPrev = 0.0;
    Scalar pvPrev = 0.0;
    Scalar parPrev = 0.0;
    Scalar pdPrev = 0.0;

    Scalar yPrevPrev = 0.0;

    Scalar yMid = 0.0;
    Scalar pvMid = 0.0;
    Scalar parMid = 0.0;
    Scalar pdMid = 0.0;
    Scalar vel = 0.0;

    Scalar sqrtC = 0.0;
    Scalar C = 0.0;
    Scalar strain1D = 0.0;
    Scalar diffGreen = 0.0;

    Scalar stressPassive = 0.0;
    Scalar diffStressPassive = 0.0;
    Scalar stressViscous = 0.0;
    Scalar diffStressViscous = 0.0;

    Scalar pAtCur = 0.0;
    Scalar pAtPrev = 0.0;
    Scalar pSvMid = 0.0;

    Scalar cavityFluxCur = 0.0;
    Scalar cavityFluxPrev = 0.0;
    Scalar dCavityFluxCur_dPv = 0.0;
    Scalar dCavityFluxCur_dPar = 0.0;

    Scalar windkesselOutflow = 0.0;
    Scalar dWindkesselOutflow_dPv = 0.0;
    Scalar dWindkesselOutflow_dPar = 0.0;

    LocalActiveDataT<Scalar> active;
  };
}

#endif
