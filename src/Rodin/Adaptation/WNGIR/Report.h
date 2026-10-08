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
      enum class Reason
      {
        IterationLimit,
        MissingInterface,
        EmptyInterface,
        DegenerateGradient,
        InvalidInitialGeometry,
        InvalidGeometry,
        GeometricTarget,
        PredictorFailure,
        NonDescentPredictor,
        InvalidScaling,
        LinearFailure,
        InnerLineSearchFailure,
        InnerIterationLimit,
        NonfiniteDirection,
        NonDescentDirection,
        LineSearchFailure,
        SmallAcceptedSteps,
        SmallEnergyChanges
      };

      /// @brief Stable log token for the typed stopping reason.
      const char* getReasonString() const
      {
        switch (reason)
        {
          case Reason::IterationLimit:
            return "iter-budget";
          case Reason::MissingInterface:
            return "missing-interface-attribute";
          case Reason::EmptyInterface:
            return "empty-interface";
          case Reason::DegenerateGradient:
            return "degenerate-target-gradient";
          case Reason::InvalidInitialGeometry:
            return "initial-state-not-strictly-feasible";
          case Reason::InvalidGeometry:
            return "geometric-validation-failed";
          case Reason::GeometricTarget:
            return "full-interface-geometric-sup-converged";
          case Reason::PredictorFailure:
            return "solve-predictor-failed";
          case Reason::NonDescentPredictor:
            return "no-descent-predictor";
          case Reason::InvalidScaling:
            return "invalid-directional-scaling";
          case Reason::LinearFailure:
            return "solve-linear-failed";
          case Reason::InnerLineSearchFailure:
            return "inner-line-search-failure";
          case Reason::InnerIterationLimit:
            return "inner-iteration-limit";
          case Reason::NonfiniteDirection:
            return "solve-nonfinite";
          case Reason::NonDescentDirection:
            return "no-descent-direction";
          case Reason::LineSearchFailure:
            return "line-search-failure";
          case Reason::SmallAcceptedSteps:
            return "best-effort-step-stagnation";
          case Reason::SmallEnergyChanges:
            return "best-effort-energy-stagnation";
        }
        return "unknown";
      }
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
      /// @brief Effective per-volume coefficient assembled for the last hinge QP.
      Real hingeCoefficient = 0;
      /// @brief Last inner Newton correction relative to the current iterate.
      Real innerRelativeCorrection = 0;
      /// @brief Euclidean norm of Mv-f+DB(v), at the last inner iterate.
      Real innerResidual = std::numeric_limits<Real>::infinity();
      /// @brief Residual divided by the fixed fitting-force norm.
      Real innerRelativeResidual = std::numeric_limits<Real>::infinity();
      Real innerResidualTolerance = 0;
      /// @brief Step factor accepted by the last inner correction.
      Real lastInnerAlpha = 0;
      /// @brief Smallest step factor accepted by the final inner solve.
      Real minInnerAlpha = 1;
      /// @brief Number of full inner Newton steps in the final inner solve.
      std::size_t fullInnerSteps = 0;
      /// @brief Number of inner Newton corrections accumulated by the solve.
      std::size_t innerIterations = 0;
      /// @brief Number of inner Newton corrections in the final outer step.
      std::size_t lastInnerIterations = 0;
      /// @brief Largest number of Newton corrections in any outer iteration.
      std::size_t maxInnerIterations = 0;
      /// @brief Whether the final inner solve met its tolerance.
      bool innerConverged = false;
      /// @brief Number of inner merit backtracks accumulated over the solve.
      std::size_t innerBacktracks = 0;
      /// @brief Maximum sampled component magnitude of the accepted physical displacement.
      Real acceptedStep = 0;
      /// @brief Minimum sampled Jacobian determinant.
      Real minJ = 1;
      /// @brief Maximum sampled Jacobian determinant.
      Real maxJ = 1;
      /// @brief Maximum sampled relative distortion.
      Real maxQRel = 1;
      /// @brief Complete-interface RMS residual.
      Real residualRMS = 0;
      /// @brief Complete-interface supremum residual.
      Real residualSup = 0;
      /// @brief RMS normalized level-set residual over the complete fitted interface.
      Real geometricRMS = std::numeric_limits<Real>::infinity();
      /// @brief Maximum sampled normalized level-set residual over the complete fitted interface.
      Real geometricSup = std::numeric_limits<Real>::infinity();
      /**
       * @brief Sampled geometric score @f$C=D_\infty/h_0^{p+1}@f$, using the
       * fixed background scale and interface FE order, not the stopping tolerance.
       * This score is not a certified Hausdorff bound or an observed error order.
       */
      Real geometricConstant = std::numeric_limits<Real>::infinity();
      /// @brief RMS unoriented normal discrepancy over the complete fitted interface.
      Real normalRMS = std::numeric_limits<Real>::infinity();

      Real geometricSupTarget = 0;
      bool geometricTargetReached = false;
      bool qualityBudgetSatisfied = false;
      /// @brief Measure of the complete interface quadrature set.
      Real interfaceMeasure = 0;
      /// @brief RMS jump of the normal field across the interface.
      Real normalJumpRMS = 0;
      /// @brief Maximum jump of the normal field across the interface.
      Real normalJumpMax = 0;
      /// @brief Final Welsch fitting energy.
      Real energy = 0;
      /// @brief Why fitting stopped; only GeometricTarget certifies the target hit.
      Reason reason = Reason::IterationLimit;
      // Wall-clock breakdown (seconds, accumulated over iterations).
      Real tAssembly = 0; ///< WNGIR variational problem assembly.
      std::size_t inactiveHingeSkips =
        0; ///< Predictor already solves the inactive-hinge model.
      std::size_t directAnalyses = 0; ///< MUMPS symbolic analyses initiated by WNGIR.
      std::size_t directFactorizations =
        0; ///< MUMPS numeric factorizations initiated by WNGIR.
      Real tSetup = 0; ///< WNGIR geometry/sigma/validation tabulation.
      Real tSolve = 0; ///< Predictor and inner linear solves.
      Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
      Real tInnerLineSearch = 0; ///< fixed-inner-merit evaluation and backtracking.
      Real tInnerAssembly = 0; ///< Inner direction-system assembly.
      Real tInnerSolve = 0; ///< Inner linear solves, excluding the predictor.
      std::size_t linearIterations = 0; ///< Accumulated linear iterations.
      std::size_t linearSolveCount = 0; ///< Number of linear solves performed.
      std::size_t maxLinearIterations = 0; ///< Largest iteration count of one solve.
      Real linearError = 0; ///< Last linear solver residual/error estimate.
  };
}

#endif
