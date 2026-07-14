/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRREPORT_H
#define RODIN_ADAPTATION_WNGIRREPORT_H

#include "WNGIRCommon.h"

namespace Rodin::Adaptation
{
  /// @brief Diagnostics produced by a WNGIR solve.
  struct WNGIRReport
  {
    /// @brief Number of nonlinear iterations performed.
      std::size_t iterations = 0;
    /// @brief Welsch scale used for the interface residual.
      Real sigma = 0;
    /// @brief Last accepted line-search factor.
      Real lastAlpha = 0;
    /// @brief Norm or magnitude of the last accepted step.
      Real acceptedStep = 0;
    /// @brief Minimum sampled Jacobian determinant.
      Real minJ = 1;
    /// @brief Maximum sampled relative distortion.
      Real maxQRel = 1;
    /// @brief Active-interface RMS residual.
      Real activeRMS = 0;
    /// @brief Active-interface supremum residual.
      Real activeSup = 0;
    /// @brief Fraction of interface quadrature in the active set.
      Real activeFraction = 0;
    /// @brief Effective RMS-over-h stopping tolerance.
      Real effectiveRMSOverHTol = 0;
    /// @brief Effective sup-over-h stopping tolerance.
      Real effectiveSupOverHTol = 0;
    /// @brief RMS jump of the normal field across the interface.
      Real normalJumpRMS = 0;
    /// @brief Maximum jump of the normal field across the interface.
      Real normalJumpMax = 0;
    /// @brief Final WNGIR energy.
      Real energy = 0;
    /// @brief Textual reason the iteration stopped.
      const char* exitReason = "iter-budget";
    // Wall-clock breakdown (seconds, accumulated over iterations).
      Real tAssembly = 0;   ///< WNGIR variational problem assembly.
      Real tSetup = 0;      ///< WNGIR geometry/sigma/validation tabulation.
      Real tBulk = 0;       ///< One-time constant bulk metric assembly.
      Real tFactor = 0;     ///< CG setup/preconditioner.
      Real tSolve = 0;      ///< CG iterations.
      Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
      std::size_t linearIterations = 0; ///< Accumulated linear iterations.
      Real linearError = 0; ///< Last linear solver residual/error estimate.
      std::size_t andersonTried = 0; ///< Number of Anderson trials.
      std::size_t andersonAccepted = 0; ///< Number of accepted Anderson trials.
      Real lastAndersonTheta = 0; ///< Last Anderson damping parameter.
  };
}

#endif
