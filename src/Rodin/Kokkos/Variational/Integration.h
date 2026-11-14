/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_INTEGRATION_H
#define RODIN_KOKKOS_VARIATIONAL_INTEGRATION_H

/**
 * @file
 * @brief Contribution storage and layout shared by Kokkos integrators.
 */
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>
#include "Rodin/Types.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Kokkos/Context/Kokkos.h"

namespace Rodin::Kokkos::Variational
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

  /**
   * @brief Arrangement of the entries of a device integrator.
   *
   * Integrators with equal layouts list the same coordinates in the same
   * order, so their values can be summed entry by entry.
   */
  struct Layout
  {
    /// @brief Integrated mesh.
    const void* mesh = nullptr;

    /// @brief Vector dimension of the trial space, zero for linear forms.
    Index trialVDim = 0;

    /// @brief Vector dimension of the test space.
    Index testVDim = 0;

    /// @brief Selected attributes; empty selects every polytope of the region.
    std::vector<::Rodin::Geometry::Attribute> attributes;

    /// @brief Integration region.
    ::Rodin::Geometry::Region region = ::Rodin::Geometry::Region::Cells;

    /// @brief Offset between the components of the trial space, zero for linear forms.
    Index trialStride = 0;

    /// @brief Offset between the components of the test space.
    Index testStride = 0;

    /// @brief Maximum scalar local count of the trial space, zero for linear forms.
    Index trialCapacity = 0;

    /// @brief Maximum scalar local count of the test space.
    Index testCapacity = 0;

    /// @brief Whether the entries couple pairs of polytopes, as global integrators do.
    bool global = false;

    /**
     * @brief Selected trial attributes of a global integrator, the
     * attributes and the region above selecting its test polytopes.
     */
    std::vector<::Rodin::Geometry::Attribute> trialAttributes = {};

    /**
     * @brief Compares two layouts member by member.
     * @param other Layout to compare with.
     * @return Whether the layouts are equal.
     */
    bool operator==(const Layout& other) const = default;
  };

  /**
   * @brief Finds where an entry of a device integrator is stored.
   * @param rank Rank of the entry, the entries of each polytope following
   * the entries of the previous polytope.
   * @param count Number of entries.
   * @param spacing Spacing of the entries, as given by the integrator.
   * @return Position of the entry past the offset of the integrator.
   *
   * When every integrated polytope owns the same number n of entries and
   * there are M > 1 polytopes, entry j of the k-th polytope is stored at
   * k + j * M and the spacing is M, so that the threads of consecutive
   * polytopes write contiguous memory. Otherwise the spacing is one and
   * the entries are stored in rank order.
   */
  KOKKOS_INLINE_FUNCTION size_t getEntryPosition(
    size_t rank, size_t count, size_t spacing)
  {
    const size_t n = count / spacing;
    return rank / n + (rank % n) * spacing;
  }

  /**
   * @brief Finds the rank of the entry stored at a position, the inverse
   * of getEntryPosition().
   * @param position Position of the entry past the offset of the integrator.
   * @param count Number of entries.
   * @param spacing Spacing of the entries, as given by the integrator.
   * @return Rank of the entry.
   */
  KOKKOS_INLINE_FUNCTION size_t getEntryRank(
    size_t position, size_t count, size_t spacing)
  {
    const size_t n = count / spacing;
    return (position % spacing) * n + position / spacing;
  }

  /**
   * @brief Coordinate list of a vector in device memory.
   *
   * Duplicate entries are summed by the consumer.
   *
   * @tparam MemorySpace Kokkos memory space of the arrays.
   */
  template <class MemorySpace>
  struct Entries
  {
    /// @brief Entry indices, test degrees of freedom.
    ::Kokkos::View<Index*, MemorySpace> indices;

    /// @brief Entry values.
    ::Kokkos::View<Real*, MemorySpace> values;
  };

}

#endif
