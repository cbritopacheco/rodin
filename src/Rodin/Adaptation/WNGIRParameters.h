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
  /// @brief Surface order including the non-polynomial composed level set.
  /// The minimum resolves the coarse curved-interface integration regression;
  /// it is not an exactness guarantee for arbitrary analytic coefficients.
  inline std::size_t wngirInterfaceQuadratureOrder(std::size_t feOrder)
  {
    return std::max<std::size_t>(12, 2 * feOrder + 2);
  }

  /// @brief Independent geometric-validation order for a finite-element order.
  inline std::size_t wngirGeometricValidationOrder(std::size_t feOrder)
  {
    return std::max<std::size_t>(14, 2 * feOrder + 4);
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
      /// Robust directional scaling of the inner model, omitting the level-set Hessian.
      bool directionalNewton = true;
      Real directionalNewtonMaxStepOverH = 1; ///< Maximum predictor motion divided by h.
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
      /// Full-interface sampled normalized-residual target; zero selects h^(p+1).
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
        Real(1e-3); ///< Relative stationarity-residual tolerance for the inner QP.
      Real primalBarrierAbsoluteTolerance = Real(1e-12); ///< Absolute inner residual tolerance.
      Real muHat = Real(90); ///< @f$\widehat\mu@f$, dimensionless
        ///< barrier/model-decrease ratio.
      Real omegaMin = 0.1; ///< @f$\omega_{\min}@f$, active-set threshold on ω.
      size_t maxBacktracks = 32; ///< Maximum halvings of a physical trial increment.
      Real armijoCoefficient =
        Real(1e-4); ///< @f$c_A@f$, Armijo sufficient-decrease coefficient.
      Real jMinRatio = 1e-8; ///< @f$j_{\min}@f$, hard inadmissibility floor.
      /// @brief Jacobian floor ratio @f$j_{\mathrm{ls}}@f$ enforced by line
      /// search.
      Real jLineSearchRatio = 1e-2;
      Real energyStagTol = 1e-8; ///< Relative energy stagnation tolerance.
      Real stepTol = 0; ///< Absolute physical accepted-displacement tolerance.
      Real acceptedStepOverHTol =
        Real(5e-4); ///< >0 stops best-effort when accepted step/h is small.
      std::size_t stagnationIterations = 5; ///< Consecutive small steps or energy changes.
      Real cgRelativeTolerance =
        1e-6; ///< @f$\tau_{\mathrm{lin}}@f$, relative residual tolerance for CG.
      std::size_t cgMaxIterations =
        1000; ///< Maximum iterations for each CG linear solve.
      std::size_t maxIterations = 30; ///< Maximum nonlinear WNGIR iterations.
      std::size_t quadratureOrder =
        0; ///< @f$p_{\mathrm{quad}}@f$ override; zero selects automatic orders.
      std::size_t geometricValidationOrder =
        0; ///< Geometric-response order; zero selects max(14, 2*(FE order) + 4).
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
