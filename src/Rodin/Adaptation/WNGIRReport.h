/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRREPORT_H
#define RODIN_ADAPTATION_WNGIRREPORT_H

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
    /// @brief Natural-gradient norm obtained from the unconstrained predictor.
      Real stationarityNorm = 0;
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
    /// @brief Last primal-barrier Newton correction relative to the predictor.
      Real primalBarrierRelativeCorrection = 0;
      /// @brief Number of primal-barrier Newton corrections accumulated by the solve.
      std::size_t primalBarrierIterations = 0;
      /// @brief Number of primal-barrier Newton corrections in the final outer step.
      std::size_t lastPrimalBarrierIterations = 0;
      /// @brief Whether the final primal-barrier inner solve met its tolerance.
      bool primalBarrierConverged = false;
      /// @brief Norm or magnitude of the last accepted step.
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
      /// @brief Measure of the active interface quadrature set.
      Real activeMeasure = 0;
      /// @brief Measure of the complete interface quadrature set.
      Real interfaceMeasure = 0;
      /// @brief Smallest observation eigenvalue on the rigid-motion space.
      Real rigidModeCoercivity = 0;
      /// @brief Smallest-to-largest rigid-mode observation eigenvalue ratio.
      Real rigidModeCoercivityRatio = 0;

      /// @brief Fraction of the accepted displacement carried by translation.
      Real rigidTranslationFraction = 0;

      /// @brief Fraction of the accepted displacement carried by rotation.
      Real rigidRotationFraction = 0;
      /// @brief Dimension of the uncontrolled rigid-motion space.
      std::size_t rigidModeDimension = 0;
      /// @brief Effective RMS-over-(h times level-set gradient) tolerance.
      Real effectiveTauRmsH = 0;
      /// @brief Effective sup-over-(h times level-set gradient) tolerance.
      Real effectiveTauInfH = 0;
      /// @brief Effective RMS tolerance in level-set units.
      Real effectiveTauRms = 0;
      /// @brief Effective supremum tolerance in level-set units.
      Real effectiveTauInf = 0;
      /// @brief RMS jump of the normal field across the interface.
      Real normalJumpRMS = 0;
      /// @brief Maximum jump of the normal field across the interface.
      Real normalJumpMax = 0;
      /// @brief Final WNGIR energy.
      Real energy = 0;
      /// @brief Textual reason the iteration stopped.
      const char* exitReason = "iter-budget";
      // Wall-clock breakdown (seconds, accumulated over iterations).
      Real tAssembly = 0; ///< WNGIR variational problem assembly.
      Real tSetup = 0; ///< WNGIR geometry/sigma/validation tabulation.
      Real tBulk = 0; ///< One-time constant bulk metric assembly.
      Real tFactor = 0; ///< CG setup/preconditioner.
      Real tSolve = 0; ///< CG iterations.
      Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
      std::size_t linearIterations = 0; ///< Accumulated linear iterations.
      Real linearError = 0; ///< Last linear solver residual/error estimate.
  };
}

#endif
