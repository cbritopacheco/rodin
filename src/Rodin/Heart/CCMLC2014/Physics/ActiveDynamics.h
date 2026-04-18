/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveDynamics.h
 * @brief Local active dynamics and linearization terms for the CCMLC2014 model.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_ACTIVEDYNAMICS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_ACTIVEDYNAMICS_H

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Piecewise-linear recruitment factor @f$ n_0(e_c) @f$ from CCMLC2014.
   *
   * @tparam Scalar Scalar numeric type.
   * @param[in] fiberDeformation Previous contractile deformation @f$ e_c @f$.
   * @returns Recruitment fraction in @f$ [0, 1] @f$.
   */
  template <class Scalar>
  inline Scalar computeRecruitmentFractionPiecewise(Scalar fiberDeformation)
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

    Scalar recruitmentFraction = Scalar(0.0);

    if (fiberDeformation < x2)
      recruitmentFraction = ((y2 - y1) / (x2 - x1)) * (fiberDeformation - x2) + y2;
    else if (fiberDeformation < x3)
      recruitmentFraction = ((y3 - y2) / (x3 - x2)) * (fiberDeformation - x3) + y3;
    else if (fiberDeformation < x4)
      recruitmentFraction = ((y4 - y3) / (x4 - x3)) * (fiberDeformation - x4) + y4;
    else if (fiberDeformation < x5)
      recruitmentFraction = y4;
    else if (fiberDeformation < x6)
      recruitmentFraction = ((y6 - y5) / (x6 - x5)) * (fiberDeformation - x6) + y6;

    return std::max<Scalar>(recruitmentFraction, Scalar(0));
  }

  /**
   * @brief Update active internal variables @f$ \gamma @f$ and @f$ \beta @f$.
   */
  template <class Input, class Scalar>
  inline void updateInternalVariables0D(
      const Input& in,
      Scalar dt,
      Scalar fiberDeformationPrevious,
      Scalar fiberDeformationCurrent,
      Scalar gammaPrevious,
      Scalar betaPrevious,
      Scalar activationDrive,
      Scalar& gammaCurrent,
      Scalar& betaCurrent,
      Scalar& recruitmentFraction)
  {
    const Scalar alpha = in.alpha;
    const Scalar k0 = in.k0;
    const Scalar sigma0 = in.sigma0;

    const Scalar activationDrivePositivePart = std::max<Scalar>(activationDrive, Scalar(0));

    recruitmentFraction = computeRecruitmentFractionPiecewise(fiberDeformationPrevious);

    const Scalar denominatorGamma =
      Scalar(1)
      + dt * (std::abs(activationDrive)
      + alpha * std::abs(fiberDeformationCurrent - fiberDeformationPrevious) / dt);

    Scalar gammaSquare =
      (gammaPrevious * gammaPrevious + dt * recruitmentFraction * k0 * activationDrivePositivePart) / denominatorGamma;
    gammaSquare = std::max<Scalar>(Scalar(1e-16), gammaSquare);
    gammaCurrent = std::sqrt(gammaSquare);

    const Scalar denominatorBeta =
      Scalar(1)
      + Scalar(0.5) * dt * recruitmentFraction * k0 * activationDrivePositivePart / gammaSquare
      + Scalar(0.5) * dt * (std::abs(activationDrive)
      + alpha * std::abs(fiberDeformationCurrent - fiberDeformationPrevious) / dt);

    betaCurrent =
      (betaPrevious + gammaCurrent * (fiberDeformationCurrent - fiberDeformationPrevious)
      + dt * recruitmentFraction * sigma0 * activationDrivePositivePart / gammaCurrent)
      / denominatorBeta;
  }

  /**
   * @brief Compute @f$\partial(\beta\gamma)/\partial e_c@f$ for local active tangent terms.
   */
  template <class Input, class Scalar>
  inline void computePartialDerivativeBetaGammaWrtFiberDeformation(
      const Input& in,
      Scalar dt,
      Scalar gammaPrevious,
      Scalar betaPrevious,
      Scalar activationDrive,
      Scalar fiberDeformationPrevious,
      Scalar fiberDeformationCurrent,
      Scalar& derivativeBetaGammaWrtFiberDeformation)
  {
    const Scalar alpha = in.alpha;
    const Scalar k0 = in.k0;
    const Scalar sigma0 = in.sigma0;

    const Scalar fiberDeformationIncrement = fiberDeformationCurrent - fiberDeformationPrevious;
    const Scalar absoluteFiberDeformationIncrement = std::abs(fiberDeformationIncrement);
    const Scalar sign =
      (fiberDeformationIncrement > Scalar(0)) ? Scalar(1) :
      (fiberDeformationIncrement < Scalar(0)) ? Scalar(-1) : Scalar(0);

    const Scalar activationDrivePositivePart = std::max<Scalar>(activationDrive, Scalar(0));
    const Scalar recruitmentFraction = computeRecruitmentFractionPiecewise(fiberDeformationPrevious);

    const Scalar Dg =
      Scalar(1) + dt * std::abs(activationDrive) + alpha * absoluteFiberDeformationIncrement;

    const Scalar Ng =
      gammaPrevious * gammaPrevious + dt * recruitmentFraction * k0 * activationDrivePositivePart;

    const Scalar gammaNewSq = std::max<Scalar>(Scalar(1e-16), Ng / Dg);
    const Scalar gammaNew = std::sqrt(gammaNewSq);

    const Scalar dGammaSq =
      -Ng * alpha * sign / (Dg * Dg);

    const Scalar dGamma =
      Scalar(0.5) * dGammaSq / gammaNew;

    const Scalar Nb =
      betaPrevious + gammaNew * fiberDeformationIncrement
      + dt * recruitmentFraction * sigma0 * activationDrivePositivePart / gammaNew;

    const Scalar Db =
      Scalar(1)
      + Scalar(0.5) * dt * recruitmentFraction * k0 * activationDrivePositivePart / gammaNewSq
      + Scalar(0.5) * dt * std::abs(activationDrive)
      + Scalar(0.5) * alpha * absoluteFiberDeformationIncrement;

    const Scalar dNb =
      dGamma * fiberDeformationIncrement
      + gammaNew
      - dt * recruitmentFraction * sigma0 * activationDrivePositivePart * dGamma / (gammaNew * gammaNew);

    const Scalar dDb =
      -Scalar(0.5) * dt * recruitmentFraction * k0 * activationDrivePositivePart
      * dGammaSq / (gammaNewSq * gammaNewSq)
      + Scalar(0.5) * alpha * sign;

    const Scalar dBeta =
      (dNb * Db - Nb * dDb) / (Db * Db);

    derivativeBetaGammaWrtFiberDeformation =
      dGamma * (Nb / Db) + gammaNew * dBeta;
  }

  /**
   * @brief Solve local active dynamics and return effective active stress/tangent.
   */
  template <class Input, class EvalData, class LocalActiveData>
  inline bool solveLocalDynamicActive(
      const Input& in,
      const EvalData& d,
      LocalActiveData& a)
  {
    using Scalar = decltype(d.y);
    a.fiberDeformationPrevious = d.sn.ec;
    a.gammaPrevious = d.sn.gamma;
    a.betaPrevious = d.sn.beta;
    a.activationDrive = in.u(d.tnp1);
    a.activationDrivePositivePart = std::max<Scalar>(a.activationDrive, Scalar(0));
    a.recruitmentFraction = computeRecruitmentFractionPiecewise(a.fiberDeformationPrevious);

    a.fiberDeformationCurrent = a.fiberDeformationPrevious;

    bool ok = false;
    for (size_t it = 0; it < in.localMaxIterations; ++it)
    {
      a.iterations = it + 1;
      a.fiberDeformationMidpoint = Scalar(0.5) * (a.fiberDeformationCurrent + a.fiberDeformationPrevious);

      updateInternalVariables0D(
          in, d.dt, a.fiberDeformationPrevious, a.fiberDeformationCurrent,
          a.gammaPrevious, a.betaPrevious, a.activationDrive,
          a.gammaCurrent, a.betaCurrent, a.recruitmentFraction);

      a.partialResidualWrtDisplacement =
        -in.Es * (Scalar(1) + Scalar(4) * d.strain1D - Scalar(2) * a.fiberDeformationMidpoint);

      {
        Scalar derivativeBetaGammaWrtFiberDeformation = 0.0;
        computePartialDerivativeBetaGammaWrtFiberDeformation(
            in, d.dt, a.gammaPrevious, a.betaPrevious, a.activationDrive,
            a.fiberDeformationPrevious, a.fiberDeformationCurrent,
            derivativeBetaGammaWrtFiberDeformation);

        a.partialResidualWrtFiberDeformation =
            Scalar(3) * std::pow(Scalar(1) + Scalar(2) * a.fiberDeformationMidpoint, 2) *
                (a.gammaCurrent * a.betaCurrent
                + in.mu * (a.fiberDeformationCurrent - a.fiberDeformationPrevious) / d.dt)
          + std::pow(Scalar(1) + Scalar(2) * a.fiberDeformationMidpoint, 3) *
                (derivativeBetaGammaWrtFiberDeformation + in.mu / d.dt)
          + Scalar(0.5) * in.Es * (Scalar(1) + Scalar(2) * d.strain1D);
      }

      const Scalar localResidual =
        (a.gammaCurrent * a.betaCurrent + in.mu * (a.fiberDeformationCurrent - a.fiberDeformationPrevious) / d.dt)
          * std::pow(Scalar(1) + Scalar(2) * a.fiberDeformationMidpoint, 3)
        - in.Es * (d.strain1D - a.fiberDeformationMidpoint) * (Scalar(1) + Scalar(2) * d.strain1D);

      a.fiberDeformationNewtonStep = localResidual / a.partialResidualWrtFiberDeformation;

      if (std::abs(localResidual) < in.localTolerance)
      {
        ok = true;
        break;
      }

      if (std::abs(a.partialResidualWrtFiberDeformation) < std::numeric_limits<Scalar>::epsilon())
        break;

      a.fiberDeformationCurrent += -in.localDamping * a.fiberDeformationNewtonStep;
    }

    a.fiberDeformationMidpoint = Scalar(0.5) * (a.fiberDeformationCurrent + a.fiberDeformationPrevious);
    updateInternalVariables0D(
        in, d.dt, a.fiberDeformationPrevious, a.fiberDeformationCurrent,
        a.gammaPrevious, a.betaPrevious, a.activationDrive,
        a.gammaCurrent, a.betaCurrent, a.recruitmentFraction);

    a.activeStressOneDimensional =
      in.Es / std::pow(Scalar(1) + Scalar(2) * a.fiberDeformationMidpoint, 2)
      * (d.strain1D - a.fiberDeformationMidpoint);

    a.partialActiveStressWrtDisplacement =
      in.Es / std::pow(Scalar(1) + Scalar(2) * a.fiberDeformationMidpoint, 2);

    a.partialActiveStressWrtFiberDeformation =
      in.Es / std::pow(Scalar(1) + Scalar(2) * a.fiberDeformationMidpoint, 3)
      * (Scalar(2) * a.fiberDeformationMidpoint - Scalar(4) * d.strain1D - Scalar(1));

    const Scalar schurComplementCorrection =
      Scalar(0.5) * a.partialActiveStressWrtFiberDeformation
      / a.partialResidualWrtFiberDeformation * a.partialResidualWrtDisplacement;

    const Scalar tangentCorrection =
      a.partialActiveStressWrtDisplacement - schurComplementCorrection;

    const Scalar rhsCorrection =
      Scalar(0.5) * a.partialActiveStressWrtFiberDeformation * a.fiberDeformationNewtonStep;

    a.activeStress = a.activeStressOneDimensional - rhsCorrection;
    a.dActiveStressWrtDisplacement = tangentCorrection * d.diffGreen;
    a.converged = ok;
    return ok;
  }
}

#endif
