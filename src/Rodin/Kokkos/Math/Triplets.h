/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MATH_TRIPLETS_H
#define RODIN_KOKKOS_MATH_TRIPLETS_H

/**
 * @file
 * @brief Sparse matrix coordinate contributions.
 */

#include <Kokkos_Core.hpp>
#include "Rodin/Types.h"

namespace Rodin::Kokkos::Math
{
  /**
   * @brief Coordinate list of a sparse matrix in device memory.
   *
   * Duplicate entries are summed by the consumer.
   *
   * @tparam MemorySpace Kokkos memory space of the arrays.
   */
  template <class MemorySpace>
  struct Triplets
  {
    /// @brief Row indices, test degrees of freedom.
    ::Kokkos::View<Index*, MemorySpace> rows;

    /// @brief Column indices, trial degrees of freedom.
    ::Kokkos::View<Index*, MemorySpace> cols;

    /// @brief Entry values.
    ::Kokkos::View<Real*, MemorySpace> values;
  };
}

#endif
