/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRPARAMETERS_H
#define RODIN_ADAPTATION_WNGIRPARAMETERS_H

#include <algorithm>
#include <cstddef>

#include "Rodin/Configure.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /// @brief Interface-assembly order for a finite-element order.
  inline std::size_t wngirInterfaceQuadratureOrder(std::size_t feOrder)
  {
    return std::max<std::size_t>(4, 2 * feOrder + 2);
  }

  /// @brief Independent geometric-validation order for a finite-element order.
  inline std::size_t wngirGeometricValidationOrder(std::size_t feOrder)
  {
    return std::max<std::size_t>(6, 2 * feOrder + 4);
  }

  /// @brief Runtime parameters controlling WNGIR assembly and iteration.
  struct WNGIRParameters
  {
      /// Affine quadratic-hinge guard widths, relative to the identity margins.
      Real qualityGuard = Real(0.1);
      /// Independent fitting, shape and distribution weights; volume terms scale with h.
      Real kappaF = 1; ///< Fitting curvature weight.
      Real kappaS = 1; ///< Shape curvature weight.
      Real kappaD = 1; ///< Distribution (current-strain regularity) weight.
      /// Robust directional Newton seed, omitting the level-set Hessian.
      bool directionalNewton = true;
      Real directionalNewtonMaxAlpha = 100;
      enum class DirectSolver
      {
        CG,
        SparseLU,
        MUMPS
      };
#ifdef RODIN_USE_MUMPS
      DirectSolver directSolver = DirectSolver::MUMPS;
#else
      DirectSolver directSolver = DirectSolver::SparseLU;
#endif
      Index directSolverThreads = 0;
      /// Physical full-interface sampled distance target; zero preserves legacy stopping.
      Real geometricSupTolerance = 0;
      Real robustScale =
        0; ///< >0 fixes the robust scale in level-set units; zero selects it automatically.
      Real h = 0; ///< reference mesh size (required).
      Real kappaJ = 1; ///< @f$\kappa_j@f$, Jacobian hinge row weight.
      Real kappaQ = 1; ///< @f$\kappa_Q@f$, relative-distortion hinge row weight.
      Real jSafe = 1e-2; ///< @f$j_{\mathrm{safe}}@f$, barrier floor on normalised j.
      Real qMax = 10; ///< @f$Q_{\max}@f$, barrier + line-search ceiling on Q.
      std::size_t primalBarrierIterations =
        15; ///< Maximum Newton corrections of the hinge-penalized QP.
      Real primalBarrierRelativeTolerance =
        Real(1e-3); ///< Relative Newton-correction tolerance for the hinge-penalized QP.
      Real muHat = Real(90); ///< @f$\widehat\mu@f$, dimensionless
        ///< barrier/model-decrease ratio.
      Real omegaMin = 0.1; ///< @f$\omega_{\min}@f$, active-set threshold on ω.
      Real alphaMin = 1e-4; ///< @f$\alpha_{\min}@f$, line-search floor.
      Real armijoCoefficient =
        Real(1e-4); ///< @f$c_A@f$, Armijo sufficient-decrease coefficient.
      Real descentFraction =
        Real(1e-4); ///< Minimum force action relative to the predictor.
      Real directionNormFactor =
        Real(10); ///< Maximum coefficient norm relative to the predictor.
      Real jMinRatio = 1e-8; ///< @f$j_{\min}@f$, hard inadmissibility floor.
      /// @brief Jacobian floor ratio @f$j_{\mathrm{ls}}@f$ enforced by line
      /// search.
      Real jLineSearchRatio = 1e-2;
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
      Real cgRelativeTolerance =
        1e-6; ///< @f$\tau_{\mathrm{lin}}@f$, relative residual tolerance for CG.
      bool cgStrictTolerance =
        false; ///< Require the requested residual without the legacy 1e-6 floor.
      std::size_t cgMaxIterations =
        1000; ///< Maximum iterations for each CG linear solve.
      std::size_t maxIterations = 30; ///< Maximum nonlinear WNGIR iterations.
      std::size_t quadratureOrder =
        0; ///< @f$p_{\mathrm{quad}}@f$ override; zero selects automatic orders.
      std::size_t geometricValidationOrder =
        0; ///< Geometric-response order; 0 ⇒ max(6, 2·(FE order) + 4).
      bool hasInterfaceAttribute = false; ///< Whether an interface marker was configured.
      Geometry::Attribute interfaceAttribute =
        0; ///< Mesh attribute identifying interface facets.
      bool trace =
        false; ///< Print inner diagnostics and validate each accepted outer geometry.
      bool traceQualityWitness =
        false; ///< Record same-point quality predictions at each trial's limiting cell.
      /// @brief Compute the rigid-observation coercivity diagnostics.
      ///
      /// The initial and final rigid-mode states are reported but never read
      /// by the solve, and they cost a generalized eigenproblem each. Set to
      /// false to skip them when the report fields are not needed.
      bool rigidDiagnostics = true;
  };
}

#endif
