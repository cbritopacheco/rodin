/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveDynamics.h
 * @brief Local active dynamics solver for the CCMLC2014 model.
 *
 * Implements the local (subcellular) active-contraction model described in
 * Caruel et al. (2014), §2.3 and §4. The internal variables @f$ \gamma @f$
 * (active stiffness) and @f$ \beta @f$ (active stress) are updated
 * implicitly, and the local fiber-deformation equilibrium is solved via a
 * Newton iteration at each global nonlinear iterate.
 *
 * The Schur-complement condensation of the local tangent onto the global
 * displacement unknown follows the approach in §4.2 of the reference.
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
   * @brief Solves the local active dynamics and computes condensed active
   *        stress/tangent for the global balance.
   *
   * At each global Newton iteration, the local fiber deformation @f$ e_c @f$
   * and internal variables @f$ \gamma, \beta @f$ are resolved iteratively.
   * The condensed active stress and its displacement derivative are then
   * returned for assembly into the global residual and Jacobian.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class ActiveDynamicsSolver
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters (active law coefficients).
       */
      explicit ActiveDynamicsSolver(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Solve the local active dynamics at the current global iterate.
       *
       * @tparam EvalData Evaluation data structure type.
       * @tparam LocalActiveData Local active data structure type.
       * @param[in] data Global evaluation data with kinematics computed.
       * @param[out] active Local active data; populated with converged internal
       *   variables and condensed active stress/tangent.
       * @returns @c true if the local Newton iteration converged.
       */
      template <class EvalData, class LocalActiveData>
      bool evaluate(const EvalData& data, LocalActiveData& active) const
      {
        using Scalar = decltype(data.y);
        active.fiberDeformationPrevious = data.sn.ec;
        active.gammaPrevious = data.sn.gamma;
        active.betaPrevious = data.sn.beta;
        active.activationDrive = m_input.u(data.tnp1);
        active.activationDrivePositivePart =
          std::max<Scalar>(active.activationDrive, Scalar(0));
        active.recruitmentFraction =
          computeRecruitmentFraction(active.fiberDeformationPrevious);

        active.fiberDeformationCurrent = active.fiberDeformationPrevious;

        bool ok = false;
        for (size_t it = 0; it < m_input.localMaxIterations; ++it)
        {
          active.iterations = it + 1;
          active.fiberDeformationMidpoint =
            Scalar(0.5) * (active.fiberDeformationCurrent
                           + active.fiberDeformationPrevious);

          updateInternalVariables(
              data.dt,
              active.fiberDeformationPrevious,
              active.fiberDeformationCurrent,
              active.gammaPrevious,
              active.betaPrevious,
              active.activationDrive,
              active.gammaCurrent,
              active.betaCurrent,
              active.recruitmentFraction);

          active.partialResidualWrtDisplacement =
            -m_input.Es * (Scalar(1) + Scalar(4) * data.strain1D
                           - Scalar(2) * active.fiberDeformationMidpoint);

          {
            Scalar derivativeBetaGammaWrtFiberDeformation = 0.0;
            computePartialDerivativeBetaGammaWrtFiberDeformation(
                data.dt,
                active.gammaPrevious,
                active.betaPrevious,
                active.activationDrive,
                active.fiberDeformationPrevious,
                active.fiberDeformationCurrent,
                derivativeBetaGammaWrtFiberDeformation);

            const Scalar ecMid2 =
              Scalar(1) + Scalar(2) * active.fiberDeformationMidpoint;

            active.partialResidualWrtFiberDeformation =
                Scalar(3) * std::pow(ecMid2, 2)
                  * (active.gammaCurrent * active.betaCurrent
                     + m_input.mu
                       * (active.fiberDeformationCurrent
                          - active.fiberDeformationPrevious) / data.dt)
              + std::pow(ecMid2, 3)
                  * (derivativeBetaGammaWrtFiberDeformation + m_input.mu / data.dt)
              + Scalar(0.5) * m_input.Es * (Scalar(1) + Scalar(2) * data.strain1D);
          }

          const Scalar ecMid3 =
            std::pow(Scalar(1) + Scalar(2) * active.fiberDeformationMidpoint, 3);
          const Scalar localResidual =
            (active.gammaCurrent * active.betaCurrent
             + m_input.mu * (active.fiberDeformationCurrent
                             - active.fiberDeformationPrevious) / data.dt)
              * ecMid3
            - m_input.Es * (data.strain1D - active.fiberDeformationMidpoint)
              * (Scalar(1) + Scalar(2) * data.strain1D);

          active.fiberDeformationNewtonStep =
            localResidual / active.partialResidualWrtFiberDeformation;

          if (std::abs(localResidual) < m_input.localTolerance)
          {
            ok = true;
            break;
          }

          if (std::abs(active.partialResidualWrtFiberDeformation)
              < std::numeric_limits<Scalar>::epsilon())
            break;

          active.fiberDeformationCurrent +=
            -m_input.localDamping * active.fiberDeformationNewtonStep;
        }

        // --- Post-convergence: compute condensed active stress ---
        active.fiberDeformationMidpoint =
          Scalar(0.5) * (active.fiberDeformationCurrent
                         + active.fiberDeformationPrevious);

        updateInternalVariables(
            data.dt,
            active.fiberDeformationPrevious,
            active.fiberDeformationCurrent,
            active.gammaPrevious,
            active.betaPrevious,
            active.activationDrive,
            active.gammaCurrent,
            active.betaCurrent,
            active.recruitmentFraction);

        const Scalar ecMid2 =
          std::pow(Scalar(1) + Scalar(2) * active.fiberDeformationMidpoint, 2);

        active.activeStressOneDimensional =
          m_input.Es / ecMid2
          * (data.strain1D - active.fiberDeformationMidpoint);

        active.partialActiveStressWrtDisplacement =
          m_input.Es / ecMid2;

        active.partialActiveStressWrtFiberDeformation =
          m_input.Es
          / (ecMid2 * (Scalar(1) + Scalar(2) * active.fiberDeformationMidpoint))
          * (Scalar(2) * active.fiberDeformationMidpoint
             - Scalar(4) * data.strain1D - Scalar(1));

        const Scalar schurComplementCorrection =
          Scalar(0.5) * active.partialActiveStressWrtFiberDeformation
          / active.partialResidualWrtFiberDeformation
          * active.partialResidualWrtDisplacement;

        const Scalar tangentCorrection =
          active.partialActiveStressWrtDisplacement - schurComplementCorrection;

        const Scalar rhsCorrection =
          Scalar(0.5) * active.partialActiveStressWrtFiberDeformation
          * active.fiberDeformationNewtonStep;

        active.activeStress = active.activeStressOneDimensional - rhsCorrection;
        active.dActiveStressWrtDisplacement = tangentCorrection * data.diffGreen;
        active.converged = ok;
        return ok;
      }

    private:
      /**
       * @brief Piecewise-linear recruitment factor @f$ n_0(e_c) @f$.
       *
       * See Caruel et al. (2014), §2.3, length-dependence of recruitment.
       *
       * @tparam Scalar Scalar numeric type.
       * @param[in] fiberDeformation Previous contractile deformation @f$ e_c @f$.
       * @returns Recruitment fraction in @f$ [0, 1] @f$.
       */
      template <class Scalar>
      static Scalar computeRecruitmentFraction(Scalar fiberDeformationPrevious)
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

        Scalar result = Scalar(0.0);

        if (fiberDeformationPrevious < x2)
          result = ((y2 - y1) / (x2 - x1)) * (fiberDeformationPrevious - x2) + y2;
        else if (fiberDeformationPrevious < x3)
          result = ((y3 - y2) / (x3 - x2)) * (fiberDeformationPrevious - x3) + y3;
        else if (fiberDeformationPrevious < x4)
          result = ((y4 - y3) / (x4 - x3)) * (fiberDeformationPrevious - x4) + y4;
        else if (fiberDeformationPrevious < x5)
          result = y4;
        else if (fiberDeformationPrevious < x6)
          result = ((y6 - y5) / (x6 - x5)) * (fiberDeformationPrevious - x6) + y6;

        return std::max<Scalar>(result, Scalar(0));
      }

      /**
       * @brief Update active internal variables @f$ \gamma @f$ and @f$ \beta @f$.
       *
       * Implements the implicit time-discrete update for the active stiffness
       * and stress internal variables (Caruel et al. 2014, §2.3).
       */
      template <class Scalar>
      void updateInternalVariables(
          Scalar dt,
          Scalar fiberDeformationPrevious,
          Scalar fiberDeformationCurrent,
          Scalar gammaPrevious,
          Scalar betaPrevious,
          Scalar activationDrive,
          Scalar& gammaCurrent,
          Scalar& betaCurrent,
          Scalar& recruitmentFraction) const
      {
        const Scalar alpha = m_input.alpha;
        const Scalar k0 = m_input.k0;
        const Scalar sigma0 = m_input.sigma0;

        const Scalar activationDrivePositivePart =
          std::max<Scalar>(activationDrive, Scalar(0));

        recruitmentFraction =
          computeRecruitmentFraction(fiberDeformationPrevious);

        const Scalar denominatorGamma =
          Scalar(1)
          + dt * (std::abs(activationDrive)
                  + alpha * std::abs(fiberDeformationCurrent
                                     - fiberDeformationPrevious) / dt);

        Scalar gammaSquare =
          (gammaPrevious * gammaPrevious
           + dt * recruitmentFraction * k0 * activationDrivePositivePart)
          / denominatorGamma;
        gammaSquare = std::max<Scalar>(Scalar(1e-16), gammaSquare);
        gammaCurrent = std::sqrt(gammaSquare);

        const Scalar denominatorBeta =
          Scalar(1)
          + Scalar(0.5) * dt * recruitmentFraction * k0
            * activationDrivePositivePart / gammaSquare
          + Scalar(0.5) * dt
            * (std::abs(activationDrive)
               + alpha * std::abs(fiberDeformationCurrent
                                  - fiberDeformationPrevious) / dt);

        betaCurrent =
          (betaPrevious
           + gammaCurrent * (fiberDeformationCurrent - fiberDeformationPrevious)
           + dt * recruitmentFraction * sigma0
             * activationDrivePositivePart / gammaCurrent)
          / denominatorBeta;
      }

      /**
       * @brief Compute @f$ \partial(\beta\gamma)/\partial e_c @f$ for
       *        the local active tangent contribution.
       *
       * This derivative is needed for the Schur-complement condensation
       * of the local fiber-deformation unknown onto the global system.
       */
      template <class Scalar>
      void computePartialDerivativeBetaGammaWrtFiberDeformation(
          Scalar dt,
          Scalar gammaPrevious,
          Scalar betaPrevious,
          Scalar activationDrive,
          Scalar fiberDeformationPrevious,
          Scalar fiberDeformationCurrent,
          Scalar& derivativeBetaGammaWrtFiberDeformation) const
      {
        const Scalar alpha = m_input.alpha;
        const Scalar k0 = m_input.k0;
        const Scalar sigma0 = m_input.sigma0;

        const Scalar fiberDeformationIncrement =
          fiberDeformationCurrent - fiberDeformationPrevious;
        const Scalar absoluteFiberDeformationIncrement =
          std::abs(fiberDeformationIncrement);
        const Scalar sign =
          (fiberDeformationIncrement > Scalar(0)) ? Scalar(1) :
          (fiberDeformationIncrement < Scalar(0)) ? Scalar(-1) : Scalar(0);

        const Scalar activationDrivePositivePart =
          std::max<Scalar>(activationDrive, Scalar(0));
        const Scalar recruitmentFraction =
          computeRecruitmentFraction(fiberDeformationPrevious);

        const Scalar Dg =
          Scalar(1) + dt * std::abs(activationDrive)
          + alpha * absoluteFiberDeformationIncrement;

        const Scalar Ng =
          gammaPrevious * gammaPrevious
          + dt * recruitmentFraction * k0 * activationDrivePositivePart;

        const Scalar gammaNewSq =
          std::max<Scalar>(Scalar(1e-16), Ng / Dg);
        const Scalar gammaNew = std::sqrt(gammaNewSq);

        const Scalar dGammaSq = -Ng * alpha * sign / (Dg * Dg);
        const Scalar dGamma = Scalar(0.5) * dGammaSq / gammaNew;

        const Scalar Nb =
          betaPrevious + gammaNew * fiberDeformationIncrement
          + dt * recruitmentFraction * sigma0
            * activationDrivePositivePart / gammaNew;

        const Scalar Db =
          Scalar(1)
          + Scalar(0.5) * dt * recruitmentFraction * k0
            * activationDrivePositivePart / gammaNewSq
          + Scalar(0.5) * dt * std::abs(activationDrive)
          + Scalar(0.5) * alpha * absoluteFiberDeformationIncrement;

        const Scalar dNb =
          dGamma * fiberDeformationIncrement
          + gammaNew
          - dt * recruitmentFraction * sigma0
            * activationDrivePositivePart * dGamma
            / (gammaNew * gammaNew);

        const Scalar dDb =
          -Scalar(0.5) * dt * recruitmentFraction * k0
            * activationDrivePositivePart
            * dGammaSq / (gammaNewSq * gammaNewSq)
          + Scalar(0.5) * alpha * sign;

        const Scalar dBeta = (dNb * Db - Nb * dDb) / (Db * Db);

        derivativeBetaGammaWrtFiberDeformation =
          dGamma * (Nb / Db) + gammaNew * dBeta;
      }

      const Input& m_input;
  };
}

#endif
