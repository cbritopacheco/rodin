/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIR_REPORT_H
#define RODIN_ADAPTATION_WNGIR_REPORT_H

#include <limits>

#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /// @brief Diagnostics produced by a WNGIR solve.
  struct WNGIRReport
  {
    /// @brief Number of nonlinear iterations performed.
      std::size_t iterations = 0;
    /// @brief Robust-loss scale used for the interface residual.
      Real sigma = 0;
    /// @brief Maximum sampled target level-set gradient on the interface.
      Real levelSetGradientScale = 0;
      /// @brief Signed fitting-force action on the unconstrained predictor; not stationarity.
      Real predictorAction = 0;
      /// @brief Directional scaling applied before constructing the hinge model.
      Real predictorScale = 1;
      /// @brief Similarity modes assigned zero by the last linear solve, not by the metric.
      size_t unresolvedSimilarityModes = 0;
      /// @brief Action of the negative energy derivative on the accepted direction.
      Real directionAction = 0;
      /// @brief Direction action divided by the unconstrained predictor action.
      Real descentRatio = 0;
      /// @brief Direction coefficient norm divided by the predictor norm.
      Real directionNormRatio = 0;
      /// @brief Actual energy decrease divided by its linear prediction.
      Real actualPredictedDecrease = 0;
      /// @brief Number of outer backtracks accumulated by the solve.
      std::size_t backtracks = 0;
      /// @brief Trials rejected by the sampled Jacobian condition.
      std::size_t jacobianRejections = 0;
      /// @brief Trials rejected by the sampled relative-distortion condition.
      std::size_t distortionRejections = 0;
      /// @brief Geometrically admissible trials rejected by the energy condition.
      std::size_t energyRejections = 0;
      /// @brief Last accepted line-search factor.
      Real lastAlpha = 0;
      /// @brief Effective per-volume coefficient assembled for the last barrier QP.
      Real primalBarrierCoefficient = 0;
      /// @brief Last primal-barrier Newton correction relative to the current iterate.
      Real primalBarrierRelativeCorrection = 0;
      /// @brief Euclidean norm of Mv-f+DB(v), at the last inner iterate.
      Real primalBarrierResidual = std::numeric_limits<Real>::infinity();
      /// @brief Residual divided by the fixed fitting-force norm.
      Real primalBarrierRelativeResidual = std::numeric_limits<Real>::infinity();
      Real primalBarrierResidualTolerance = 0;
      /// @brief Step factor accepted by the last primal-barrier correction.
      Real lastPrimalBarrierAlpha = 0;
      /// @brief Smallest step factor accepted by the final primal-barrier solve.
      Real minPrimalBarrierAlpha = 1;
      /// @brief Number of full primal-barrier Newton steps in the final inner solve.
      std::size_t fullPrimalBarrierSteps = 0;
      /// @brief Number of primal-barrier Newton corrections accumulated by the solve.
      std::size_t primalBarrierIterations = 0;
      /// @brief Number of primal-barrier Newton corrections in the final outer step.
      std::size_t lastPrimalBarrierIterations = 0;
      /// @brief Largest number of Newton corrections in any outer iteration.
      std::size_t maxPrimalBarrierIterations = 0;
      /// @brief Whether the final primal-barrier inner solve met its tolerance.
      bool primalBarrierConverged = false;
      /// @brief Number of inner merit backtracks accumulated over the solve.
      std::size_t primalBarrierBacktracks = 0;
      /// @brief Maximum sampled component magnitude of the accepted physical displacement.
      Real acceptedStep = 0;
      /// @brief Minimum sampled Jacobian determinant.
      Real minJ = 1;
      /// @brief Maximum sampled Jacobian determinant.
      Real maxJ = 1;
      /// @brief Maximum sampled relative distortion.
      Real maxQRel = 1;
      /// @brief Active-interface RMS residual.
      Real activeRMS = 0;
      /// @brief Active-interface supremum residual.
      Real activeSup = 0;
      /// @brief Fraction of interface quadrature in the active set.
      Real activeFraction = 0;
      /// @brief RMS normalized level-set residual over the complete fitted interface.
      Real geometricRMS = std::numeric_limits<Real>::infinity();
      /// @brief Maximum normalized level-set residual over the complete fitted interface.
      Real geometricSup = std::numeric_limits<Real>::infinity();
      /// @brief RMS unoriented normal discrepancy over the complete fitted interface.
      Real normalRMS = std::numeric_limits<Real>::infinity();

      Real geometricSupTarget = 0;
      bool geometricTargetReached = false;
      bool qualityBudgetSatisfied = false;
      /// @brief Measure of the active interface quadrature set.
      Real activeMeasure = 0;
      /// @brief Measure of the complete interface quadrature set.
      Real interfaceMeasure = 0;
      /// @brief Smallest observation eigenvalue on the rigid-motion space.
      Real rigidModeCoercivity = 0;
      /// @brief Smallest-to-largest rigid-mode observation eigenvalue ratio.
      Real rigidModeCoercivityRatio = 0;

      /// @brief Dimension of the uncontrolled rigid-motion space.
      std::size_t rigidModeDimension = 0;
      /// @brief RMS jump of the normal field across the interface.
      Real normalJumpRMS = 0;
      /// @brief Maximum jump of the normal field across the interface.
      Real normalJumpMax = 0;
      /// @brief Final Welsch fitting energy.
      Real energy = 0;
      /// @brief Textual reason the iteration stopped.
      const char* exitReason = "iter-budget";
      // Wall-clock breakdown (seconds, accumulated over iterations).
      Real tAssembly = 0; ///< WNGIR variational problem assembly.
      std::size_t inactiveHingeSkips =
        0; ///< Predictor already solves the inactive-hinge model.
      std::size_t directAnalyses = 0; ///< MUMPS symbolic analyses initiated by WNGIR.
      std::size_t directFactorizations =
        0; ///< MUMPS numeric factorizations initiated by WNGIR.
      Real tSetup = 0; ///< WNGIR geometry/sigma/validation tabulation.
      Real tFactor = 0; ///< Linear solver setup/preconditioner.
      Real tSolve = 0; ///< Predictor and inner linear solves.
      Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
      Real tPrimalBarrierLineSearch =
        0; ///< fixed-inner-merit evaluation and backtracking.
      Real tPrimalBarrierAssembly = 0; ///< Inner direction-system assembly.
      Real tPrimalBarrierSolve = 0; ///< Inner linear solves, excluding the predictor.
      std::size_t linearIterations = 0; ///< Accumulated linear iterations.
      std::size_t linearSolveCount = 0; ///< Number of linear solves performed.
      std::size_t maxLinearIterations = 0; ///< Largest iteration count of one solve.
      Real linearError = 0; ///< Last linear solver residual/error estimate.
  };
}

#endif
