/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_ACTIVEDYNAMICS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_ACTIVEDYNAMICS_H

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>

namespace Rodin::Heart::CCMLC2014::Physics
{
  template <class Scalar>
  inline Scalar n0_piecewise(Scalar fib0)
  {
    const Scalar x1 = Scalar(-0.4);
    const Scalar y1 = Scalar(0.0);
    const Scalar x2 = Scalar(0.3);
    const Scalar y2 = Scalar(0.38);
    const Scalar x3 = Scalar(0.73);
    const Scalar y3 = Scalar(0.74);
    const Scalar x4 = Scalar(1.0);
    const Scalar y4 = Scalar(1.0);
    const Scalar x5 = Scalar(1.3);
    const Scalar y5 = Scalar(1.0);
    const Scalar x6 = Scalar(2.4);
    const Scalar y6 = Scalar(0.0);

    Scalar n0 = Scalar(0.0);

    if (fib0 < x2)
      n0 = ((y2 - y1) / (x2 - x1)) * (fib0 - x2) + y2;
    else if (fib0 < x3)
      n0 = ((y3 - y2) / (x3 - x2)) * (fib0 - x3) + y3;
    else if (fib0 < x4)
      n0 = ((y4 - y3) / (x4 - x3)) * (fib0 - x4) + y4;
    else if (fib0 < x5)
      n0 = y4;
    else if (fib0 < x6)
      n0 = ((y6 - y5) / (x6 - x5)) * (fib0 - x6) + y6;

    return std::max<Scalar>(n0, Scalar(0));
  }

  template <class Input, class Scalar>
  inline void updateInternalVariables0D(
      const Input& in,
      Scalar dt,
      Scalar fib0,
      Scalar fib1,
      Scalar gammaOld,
      Scalar betaOld,
      Scalar u1,
      Scalar& gammaNew,
      Scalar& betaNew,
      Scalar& n0)
  {
    const Scalar alpha = in.alpha;
    const Scalar k0 = in.k0;
    const Scalar sigma0 = in.sigma0;

    const Scalar u1Plus = std::max<Scalar>(u1, Scalar(0));

    n0 = n0_piecewise(fib0);

    const Scalar denominatorGamma =
      Scalar(1) + dt * (std::abs(u1) + alpha * std::abs(fib1 - fib0) / dt);

    Scalar gammaSquare =
      (gammaOld * gammaOld + dt * n0 * k0 * u1Plus) / denominatorGamma;
    gammaSquare = std::max<Scalar>(Scalar(1e-16), gammaSquare);
    gammaNew = std::sqrt(gammaSquare);

    const Scalar denominatorBeta =
      Scalar(1)
      + Scalar(0.5) * dt * n0 * k0 * u1Plus / gammaSquare
      + Scalar(0.5) * dt * (std::abs(u1) + alpha * std::abs(fib1 - fib0) / dt);

    betaNew =
      (betaOld + gammaNew * (fib1 - fib0) + dt * n0 * sigma0 * u1Plus / gammaNew)
      / denominatorBeta;
  }

  template <class Input, class Scalar>
  inline void getPartialDerivativesInternalVariablesWrtFibDeformation(
      const Input& in,
      Scalar dt,
      Scalar gammaOld,
      Scalar betaOld,
      Scalar u1,
      Scalar fib0,
      Scalar fib1,
      Scalar& derBetaGammaWrtFib)
  {
    const Scalar alpha = in.alpha;
    const Scalar k0 = in.k0;
    const Scalar sigma0 = in.sigma0;

    const Scalar deltaFib = fib1 - fib0;
    const Scalar absDelta = std::abs(deltaFib);
    const Scalar sign =
      (deltaFib > Scalar(0)) ? Scalar(1) :
      (deltaFib < Scalar(0)) ? Scalar(-1) : Scalar(0);

    const Scalar u1Plus = std::max<Scalar>(u1, Scalar(0));
    const Scalar n0 = n0_piecewise(fib0);

    const Scalar Dg =
      Scalar(1) + dt * std::abs(u1) + alpha * absDelta;

    const Scalar Ng =
      gammaOld * gammaOld + dt * n0 * k0 * u1Plus;

    const Scalar gammaNewSq = std::max<Scalar>(Scalar(1e-16), Ng / Dg);
    const Scalar gammaNew = std::sqrt(gammaNewSq);

    const Scalar dGammaSq =
      -Ng * alpha * sign / (Dg * Dg);

    const Scalar dGamma =
      Scalar(0.5) * dGammaSq / gammaNew;

    const Scalar Nb =
      betaOld + gammaNew * deltaFib + dt * n0 * sigma0 * u1Plus / gammaNew;

    const Scalar Db =
      Scalar(1)
      + Scalar(0.5) * dt * n0 * k0 * u1Plus / gammaNewSq
      + Scalar(0.5) * dt * std::abs(u1)
      + Scalar(0.5) * alpha * absDelta;

    const Scalar dNb =
      dGamma * deltaFib
      + gammaNew
      - dt * n0 * sigma0 * u1Plus * dGamma / (gammaNew * gammaNew);

    const Scalar dDb =
      -Scalar(0.5) * dt * n0 * k0 * u1Plus * dGammaSq / (gammaNewSq * gammaNewSq)
      + Scalar(0.5) * alpha * sign;

    const Scalar dBeta =
      (dNb * Db - Nb * dDb) / (Db * Db);

    derBetaGammaWrtFib =
      dGamma * (Nb / Db) + gammaNew * dBeta;
  }

  template <class Input, class EvalData, class LocalActiveData>
  inline bool solveLocalDynamicActive(
      const Input& in,
      const EvalData& d,
      LocalActiveData& a)
  {
    using Scalar = decltype(d.y);
    a.fib0 = d.sn.ec;
    a.gammaOld = d.sn.gamma;
    a.betaOld = d.sn.beta;
    a.u1 = in.u(d.tnp1);
    a.u1Plus = std::max<Scalar>(a.u1, Scalar(0));
    a.n0 = n0_piecewise(a.fib0);

    a.fib1 = a.fib0;

    bool ok = false;
    for (size_t it = 0; it < in.localMaxIterations; ++it)
    {
      a.iterations = it + 1;
      a.fib12 = Scalar(0.5) * (a.fib1 + a.fib0);

      updateInternalVariables0D(
          in, d.dt, a.fib0, a.fib1, a.gammaOld, a.betaOld, a.u1,
          a.gammaNew, a.betaNew, a.n0);

      a.k21 = -in.Es * (Scalar(1) + Scalar(4) * d.strain1D - Scalar(2) * a.fib12);

      {
        Scalar derBetaGammaWrtFib = 0.0;
        getPartialDerivativesInternalVariablesWrtFibDeformation(
            in, d.dt, a.gammaOld, a.betaOld, a.u1, a.fib0, a.fib1, derBetaGammaWrtFib);

        a.k22 =
            Scalar(3) * std::pow(Scalar(1) + Scalar(2) * a.fib12, 2) *
                (a.gammaNew * a.betaNew + in.mu * (a.fib1 - a.fib0) / d.dt)
          + std::pow(Scalar(1) + Scalar(2) * a.fib12, 3) *
                (derBetaGammaWrtFib + in.mu / d.dt)
          + Scalar(0.5) * in.Es * (Scalar(1) + Scalar(2) * d.strain1D);
      }

      const Scalar Rraw =
        (a.gammaNew * a.betaNew + in.mu * (a.fib1 - a.fib0) / d.dt)
          * std::pow(Scalar(1) + Scalar(2) * a.fib12, 3)
        - in.Es * (d.strain1D - a.fib12) * (Scalar(1) + Scalar(2) * d.strain1D);

      a.krc_k22 = Rraw / a.k22;

      if (std::abs(Rraw) < in.localTolerance)
      {
        ok = true;
        break;
      }

      if (std::abs(a.k22) < std::numeric_limits<Scalar>::epsilon())
        break;

      a.fib1 += -in.localDamping * a.krc_k22;
    }

    a.fib12 = Scalar(0.5) * (a.fib1 + a.fib0);
    updateInternalVariables0D(
        in, d.dt, a.fib0, a.fib1, a.gammaOld, a.betaOld, a.u1,
        a.gammaNew, a.betaNew, a.n0);

    a.sigma1d =
      in.Es / std::pow(Scalar(1) + Scalar(2) * a.fib12, 2)
      * (d.strain1D - a.fib12);

    a.partialSigma1dWrtDisp =
      in.Es / std::pow(Scalar(1) + Scalar(2) * a.fib12, 2);

    a.partialSigma1dWrtEc =
      in.Es / std::pow(Scalar(1) + Scalar(2) * a.fib12, 3)
      * (Scalar(2) * a.fib12 - Scalar(4) * d.strain1D - Scalar(1));

    const Scalar coefSchurD2W =
      Scalar(0.5) * a.partialSigma1dWrtEc / a.k22 * a.k21;

    const Scalar tangentCorrection =
      a.partialSigma1dWrtDisp - coefSchurD2W;

    const Scalar rhsCorrection =
      Scalar(0.5) * a.partialSigma1dWrtEc * a.krc_k22;

    a.stressActive = a.sigma1d - rhsCorrection;
    a.diffStressActive = tangentCorrection * d.diffGreen;
    a.converged = ok;
    return ok;
  }
}

#endif
