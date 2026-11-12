/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_ROCSOLVERRFPIVOTS_H
#define RODIN_KOKKOS_SOLVER_ROCSOLVERRFPIVOTS_H

/**
 * @file
 * @brief Backend reduction of exact finite, nonzero factor pivot checks.
 */

#include <limits>
#include <Kokkos_Core.hpp>

#include "Rodin/Kokkos/Utility.h"

namespace Rodin::Solver::Internal
{
  /**
   * @brief Certifies that all factor pivots are finite and nonzero.
   *
   * ## Architecture
   * A Kokkos reduction counts invalid gathered pivots in their memory space.
   * Only that scalar count returns to the host. The same kernel runs with a
   * host execution space in tests and with HIP in rocSOLVER. No tolerance or
   * norm is used, so every finite nonzero pivot is accepted.
   */
  struct RocSolverRFPivots
  {
    /**
     * @brief Tests a pivot without introducing a numerical tolerance.
     * @param value Pivot value.
     * @return Whether the pivot is finite and nonzero.
     */
    KOKKOS_INLINE_FUNCTION static bool valid(double value)
    {
      constexpr double Maximum = std::numeric_limits<double>::max();
      return value != 0 && value >= -Maximum && value <= Maximum;
    }

    /**
     * @brief Checks every pivot on the supplied execution space.
     * @tparam Execution Execution space accessing the pivot view.
     * @tparam View One-dimensional view of factor pivots.
     * @param execution Execution instance on which the pivots are ready.
     * @param values Gathered pivots, which may be empty.
     * @return Whether every pivot is finite and nonzero, after reduction completion.
     */
    template <class Execution, class View>
    static bool check(const Execution& execution, const View& values)
    {
      static_assert(::Kokkos::SpaceAccessibility<Execution,
        typename View::memory_space>::accessible);
      size_t invalid = 0;
      ::Kokkos::parallel_reduce("Rodin rocSOLVER check pivots",
        ::Kokkos::RangePolicy<Execution, ::Kokkos::IndexType<size_t>>(
          execution, 0, values.extent(0)),
        Check<View>{values}, invalid);
      ::Rodin::Kokkos::Internal::checkLaunches("Rodin rocSOLVER check pivots");
      return invalid == 0;
    }

    /// @brief Counts invalid pivots without copying their values to the host.
    template <class View>
    struct Check
    {
      /// @brief Device pivots to inspect.
      View values;

      /**
       * @param i Pivot index.
       * @param invalid Number of invalid pivots in the reduction.
       */
      KOKKOS_INLINE_FUNCTION void operator()(size_t i, size_t& invalid) const
      {
        if (!valid(values(i)))
          ++invalid;
      }
    };
  };
}

#endif
