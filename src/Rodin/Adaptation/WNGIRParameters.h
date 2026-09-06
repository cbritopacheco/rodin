/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRPARAMETERS_H
#define RODIN_ADAPTATION_WNGIRPARAMETERS_H

#include <cstddef>

#include "Rodin/Geometry/Types.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /// @brief Runtime parameters controlling WNGIR assembly and iteration.
  struct WNGIRParameters
  {
      Real robustScale =
        0; ///< >0 fixes the robust scale in level-set units; zero selects it automatically.
      Real h = 0; ///< reference mesh size (required).
      Real kappaBulk =
        Real(0.001); ///< @f$\kappa_{\mathrm{bulk}}@f$, dimensionless
                           ///< bulk-strain coefficient.
      Real rDiv = 1; ///< @f$r_D@f$, divergence/deviatoric bulk-coefficient ratio.
      Real kappaObs =
        1; ///< @f$\kappa_{\mathrm{obs}}@f$, surface observation metric weight.
      Real kappaJ = 1; ///< @f$\kappa_j@f$, Jacobian barrier row weight.
      Real kappaQ = 1; ///< @f$\kappa_Q@f$, relative-distortion barrier row weight.
      Real jSafe = 1e-2; ///< @f$j_{\mathrm{safe}}@f$, barrier floor on normalised j.
      Real qMax = 10; ///< @f$Q_{\max}@f$, barrier + line-search ceiling on Q.
      std::size_t primalBarrierIterations =
        8; ///< Maximum Newton corrections of the QP barrier.
      Real primalBarrierRelativeTolerance =
        Real(1e-2); ///< Relative Newton-correction tolerance for the QP barrier.
      Real muHat = Real(0.3); ///< @f$\widehat\mu@f$, dimensionless
                              ///< barrier/model-decrease ratio.
      Real thetaBoundary = Real(0.95); ///< @f$\tau@f$, strict-feasibility fraction.
      Real omegaMin = 0.1; ///< @f$\omega_{\min}@f$, active-set threshold on ω.
      Real alphaMin = 1e-4; ///< @f$\alpha_{\min}@f$, line-search floor.
      Real armijoCoefficient =
        Real(1e-4); ///< @f$c_A@f$, Armijo sufficient-decrease coefficient.
      Real descentFraction =
        Real(1e-4); ///< Minimum force action relative to the predictor.
      Real directionNormFactor =
        Real(10); ///< Maximum coefficient norm relative to the predictor.
      Real jMinRatio = 1e-8; ///< @f$j_{\min}@f$, hard inadmissibility floor.
      Real jLineSearchRatio =
        1e-2; ///< @f$j_{\mathrm{ls}}@f$, Jacobian floor ratio enforced by
              ///< line search.
      /// @brief Absolute RMS tolerance @f$\tau_{\mathrm{rms}}@f$; zero selects
      /// four times the level-set mesh scale.
      Real tauRms = 1e-12;

      /// @brief Absolute supremum tolerance @f$\tau_\infty@f$; zero selects ten
      /// times the level-set mesh scale.
      Real tauInf = 1e-12;

      /// @brief Lower bound on the scale-aware tolerance
      /// @f$\tau_{\mathrm{rms},h}@f$.
      Real tauRmsHFloor = Real(0.005);

      /// @brief Lower bound on the scale-aware tolerance @f$\tau_{\infty,h}@f$.
      Real tauInfHFloor = 0;
      /// @brief Factor @f$\tau^{\mathrm{rms}}_{\mathrm{jump}}@f$ multiplying the
      /// normal-jump estimate in @f$\tau_{\mathrm{rms},h}@f$.
      Real tauJumpRms = 0;

      /// @brief Factor @f$\tau^{\infty}_{\mathrm{jump}}@f$ multiplying the
      /// normal-jump estimate in @f$\tau_{\infty,h}@f$.
      Real tauJumpInf = 0;
      Real energyStagTol = 1e-8; ///< Relative energy stagnation tolerance.
      Real stepTol = 0; ///< ≤0 ⇒ 1e-4·h.
      Real acceptedStepOverHTol =
        Real(5e-4); ///< >0 stops best-effort when accepted step/h is small.
      Real rigidStabilisationLevel =
        Real(1); ///< @f$\rho@f$, lifts weakly observed rigid modes to this
                   ///< fraction of the stiffest rigid mode; zero disables the
                   ///< stabilisation.
      Real cgRelativeTolerance =
        1e-6; ///< @f$\tau_{\mathrm{lin}}@f$, relative residual tolerance for CG.
      std::size_t cgMaxIterations =
        0; ///< 0 ⇒ min(2000, max(100, 2*ndofs)).
      std::size_t maxIterations = 200; ///< Maximum nonlinear WNGIR iterations.
      std::size_t quadratureOrder = 0; ///< @f$p_{\mathrm{quad}}@f$; 0 ⇒ 2·(FE order).
      bool hasInterfaceAttribute = false; ///< Whether an interface marker was configured.
      Geometry::Attribute interfaceAttribute =
        0; ///< Mesh attribute identifying interface facets.
      bool trace = false; ///< Print per-iteration diagnostics when true.
  };
}

#endif
