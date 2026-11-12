/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_ROCSOLVERRFPIVOTS_H
#define RODIN_SOLVER_ROCSOLVERRFPIVOTS_H

/**
 * @file
 * @brief Exact finite, nonzero pivot checks for rocSOLVER refactorization.
 */

#include <cstddef>
#include <limits>

#include "Rodin/Configure.h"

#ifdef RODIN_USE_ROCSOLVER
#include <hip/hip_runtime_api.h>
#endif

namespace Rodin::Solver::Internal
{
  /**
   * @brief Certifies that all factor pivots are finite and nonzero.
   *
   * ## Architecture
   * The HIP kernel checks each gathered pivot and atomically combines failures
   * into one retained device flag. Only that flag returns to the host. The
   * scalar predicate has the same exact semantics on host and device; it
   * accepts arbitrarily small nonzero pivots and does not sum their magnitudes.
   */
  struct RocSolverRFPivots
  {
    /**
     * @brief Tests a pivot without introducing a numerical tolerance.
     * @param value Pivot value.
     * @return Whether the pivot is finite and nonzero.
     */
#if defined(__HIPCC__)
    __host__ __device__
#endif
      static constexpr bool
      valid(double value)
    {
      constexpr double Maximum = std::numeric_limits<double>::max();
      return value != 0 && value >= -Maximum && value <= Maximum;
    }

#ifdef RODIN_USE_ROCSOLVER
    /**
     * @brief Checks a device array on the specified HIP stream.
     * @param values Gathered device pivots; may be null when count is zero.
     * @param count Number of pivots.
     * @param flag Retained device storage for one integer mismatch flag.
     * @param result Receives the result after successful stream completion.
     * @param stream Stream on which the gathered pivots are ready.
     * @return HIP status; a failed launch or copy never certifies the pivots.
     */
    static hipError_t check(
      const double* values, size_t count, int* flag, bool& result, hipStream_t stream);
#endif
  };
}

#endif
