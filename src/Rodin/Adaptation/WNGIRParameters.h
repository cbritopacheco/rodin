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
      Real kappaBulk = Real(0.00703125); ///< Dimensionless bulk-strain coefficient.
      Real rDiv = 1; ///< Divergence/deviatoric bulk-coefficient ratio.
      Real kappaObs = 1; ///< surface observation metric weight.
      Real tauTan =
        Real(0.05); ///< Tangential weight of the hybrid rank-one metric.
      Real kappaInit = 0; ///< Normal-offset initializer; zero gives a cold start.
      Real initialGuessCapH =
        2; ///< Cap normal-offset initialization by this multiple of h.
      Real kappaJ = 1; ///< Jacobian barrier row weight.
      Real kappaQ = 1; ///< Relative-distortion barrier row weight.
      Real jSafe = 1e-2; ///< barrier floor on j.
      Real qMax = 10; ///< barrier + line-search ceiling on Q.
      std::size_t primalBarrierIterations =
        8; ///< Maximum Newton corrections of the QP barrier.
      Real primalBarrierRelativeTolerance =
        Real(1e-2); ///< Relative Newton-correction tolerance for the QP barrier.
      bool requirePrimalBarrierConvergence =
        true; ///< Reject a primal-barrier direction without an inner certificate.
      Real muHat = Real(0.3); ///< Dimensionless barrier/model-decrease ratio.
      Real thetaBoundary = Real(0.95); ///< Strict-feasibility fraction.
      Real omegaMin = 0.1; ///< active-set threshold on ω.
      Real alphaMin = 1e-4; ///< line-search floor.
      bool admissibilityChecks = true; ///< Enforce true-geometry j and Q bounds.
      bool energyLineSearch = true; ///< Require WNGIR energy decrease in line search.
      Real armijoCoefficient = Real(1e-4); ///< Armijo sufficient-decrease coefficient.
      Real descentFraction =
        Real(1e-4); ///< Minimum force action relative to the predictor.
      Real directionNormFactor =
        Real(10); ///< Maximum coefficient norm relative to the predictor.
      Real jMinRatio = 1e-8; ///< hard inadmissibility floor.
      Real jLineSearchRatio = 1e-2; ///< Jacobian floor ratio enforced by line search.
      Real activeRMSTol = 0; ///< ≤0 ⇒ 4h² times the level-set gradient scale.
      Real activeSupTol = 0; ///< ≤0 ⇒ 10h² times the level-set gradient scale.
      Real activeRMSOverHTol =
        0; ///< >0 bounds RMS divided by h times the level-set gradient scale.
      Real activeSupOverHTol =
        0; ///< >0 bounds sup divided by h times the level-set gradient scale.
      bool geometryAwareTolerances = true; ///< Enable dimension-aware residual floors.
      Real rmsFloor2D = Real(0.05); ///< Minimum RMS residual floor in 2D, divided by h.
      Real supFloor2D = Real(0.25); ///< Minimum sup residual floor in 2D, divided by h.
      Real rmsFloor3D = Real(0.03); ///< Minimum RMS residual floor in 3D, divided by h.
      Real supFloor3D = Real(0.20); ///< Minimum sup residual floor in 3D, divided by h.
      Real rmsNormalJumpFactor =
        Real(0.03); ///< RMS tolerance factor multiplying the normal-jump estimate.
      Real supNormalJumpFactor =
        Real(0.05); ///< Sup tolerance factor multiplying the normal-jump estimate.
      Real energyStagTol = 1e-4; ///< Relative energy stagnation tolerance.
      Real stationarityTolerance = 0; ///< >0 enables natural-gradient stopping.
      Real stepTol = 0; ///< ≤0 ⇒ 1e-4·h.
      Real acceptedStepOverHTol =
        Real(5e-3); ///< >0 stops best-effort when accepted step/h is small.
      Real cgRelativeTolerance = 1e-6; ///< relative residual tolerance for CG.
      std::size_t cgMaxIterations = 0; ///< 0 ⇒ min(2000, max(100, 2*ndofs)).
      std::size_t maxIterations = 200; ///< Maximum nonlinear WNGIR iterations.
      std::size_t quadratureOrder = 0; ///< 0 ⇒ 2·(FE order).
      bool hasInterfaceAttribute = false; ///< Whether an interface marker was configured.
      Geometry::Attribute interfaceAttribute =
        0; ///< Mesh attribute identifying interface facets.
      FlatSet<Geometry::Attribute> dirichletAttributes; ///< Zero-displacement boundaries.
      bool trace = false; ///< Print per-iteration diagnostics when true.
  };
}

#endif
