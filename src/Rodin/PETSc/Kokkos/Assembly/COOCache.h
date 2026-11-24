/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_ASSEMBLY_COOCACHE_H
#define RODIN_PETSC_KOKKOS_ASSEMBLY_COOCACHE_H

/**
 * @file
 * @brief Reusable device contributions and symbolic state of PETSc COO assembly.
 */

#include "Distribution.h"
#include <petscvec.h>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>
#include "Rodin/Kokkos/Variational/Integration.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Block of the coordinate list shared by integrators with the same
   * entries.
   */
  struct COOSegment
  {
    /// @brief Arrangement of the entries.
    ::Rodin::Kokkos::Variational::Layout layout;

    /// @brief Offset of the test space in the system.
    Index row = 0;

    /// @brief Offset of the trial space in the system.
    Index col = 0;

    /// @brief Number of entries.
    size_t count = 0;

    /// @brief Spacing of the entries, as given to getEntryPosition().
    size_t spacing = 1;

    /// @brief Whether the entries are the diagonals of the component blocks only.
    bool diagonal = false;

    /**
     * @brief Compares all the members of the segments.
     * @param other Segment compared.
     * @return Whether the segments have the same entries.
     */
    bool operator==(const COOSegment& other) const = default;
  };

  /**
   * @brief Permutation sorting a coordinate list, as indices in the list.
   *
   * The indices are stored on 32 bits when the list is short enough, which
   * halves the traffic of the permutation in the sums of the values per
   * coordinate, and on the width of Index otherwise. One of the two views
   * is empty.
   *
   * @tparam MemorySpace Kokkos memory space of the indices.
   */
  template <class MemorySpace>
  struct COOOrder
  {
    /// @brief Indices, when the list has at most 2^32 - 1 entries.
    ::Kokkos::View<std::uint32_t*, MemorySpace> narrow;

    /// @brief Indices, when the list is longer.
    ::Kokkos::View<Index*, MemorySpace> wide;

    /**
     * @brief Calls the function with the view that holds the indices.
     * @param function Callable taking the view of the indices.
     * @return Result of the call.
     */
    template <class Function>
    decltype(auto) visit(Function&& function) const
    {
      if (wide.extent(0) > 0)
        return function(wide);
      return function(narrow);
    }
  };

  struct COODistribution;

  /**
   * @brief Device buffers and sparsity pattern kept between the assemblies
   * of a PETSc system.
   *
   * PETSc analyses the coordinate list in MatSetPreallocationCOO, which is
   * much more expensive than MatSetValuesCOO. The pattern is kept to skip
   * that analysis, and the coordinates are not even listed again, when the
   * same matrix is assembled again with the same segments. PETSc is given
   * the distinct coordinates sorted by row then column, which it analyses
   * in linear time, and the values summed per coordinate on the device. Integrators
   * restricted to attributes list their coordinates at every assembly,
   * since the attributes of the cells may change, and the pattern is
   * analysed again only when they moved. The buffers are reused
   * while their size does not change.
   *
   * A copy starts empty, so copies never share buffers.
   *
   * @tparam MemorySpace Kokkos memory space of the buffers.
   */
  template <class MemorySpace>
  struct COOCache
  {
    COOCache() = default;

    /**
     * @brief Constructs an empty cache, without sharing the buffers.
     * @param other Cache copied, whose buffers are not copied.
     */
    COOCache([[maybe_unused]] const COOCache& other) {}

    /**
     * @brief Empties the cache, without sharing the buffers.
     * @param other Cache assigned, whose buffers are not copied.
     * @return Reference to this cache.
     */
    COOCache& operator=(const COOCache& other)
    {
      if (this != &other)
        *this = COOCache();
      return *this;
    }

    /**
     * @brief Takes the buffers of the moved cache.
     * @param other Cache moved.
     */
    COOCache(COOCache&& other) = default;

    /**
     * @brief Takes the buffers of the moved cache.
     * @param other Cache moved.
     * @return Reference to this cache.
     */
    COOCache& operator=(COOCache&& other) = default;

    /// @brief PETSc identifier of the matrix preallocated with the pattern, or -1.
    PetscInt64 matrix = -1;

    /**
     * @brief PETSc identifier of the distributed right-hand side
     * preallocated with the local coordinates, or -1.
     */
    PetscInt64 vector = -1;

    /// @brief Whether the coordinates hold the pattern of the segments.
    bool indexed = false;

    /// @brief Segments of the pattern.
    std::vector<COOSegment> segments;

    /**
     * @brief Row coordinates, empty while they are released.
     *
     * They are kept in the PETSc index type, which receives them and takes
     * half of the listed indices when it is 32 bits wide.
     */
    ::Kokkos::View<PetscInt*, MemorySpace> rows;

    /// @brief Column coordinates, empty while they are released.
    ::Kokkos::View<PetscInt*, MemorySpace> cols;

    /// @brief Number of coordinates of the pattern.
    size_t listed = 0;

    /**
     * @brief Whether the coordinates were released, the distinct ones being
     * kept to list them again.
     */
    bool released = false;

    /// @brief Rows of the distinct coordinates, kept while they are released.
    ::Kokkos::View<PetscInt*, MemorySpace> distinctRows;

    /// @brief Columns of the distinct coordinates, kept as distinctRows.
    ::Kokkos::View<PetscInt*, MemorySpace> distinctCols;

    /// @brief Row coordinates of the restricted segments, packed as nextRows.
    ::Kokkos::View<Index*, MemorySpace> restrictedRows;

    /// @brief Column coordinates of the restricted segments, as restrictedRows.
    ::Kokkos::View<Index*, MemorySpace> restrictedCols;

    /**
     * @brief Row coordinates being listed, all of them or only those of the
     * restricted segments, packed, when they alone are listed again.
     */
    ::Kokkos::View<Index*, MemorySpace> nextRows;

    /// @brief Column coordinates being listed, as nextRows.
    ::Kokkos::View<Index*, MemorySpace> nextCols;

    /// @brief Values of the restricted segments listed again, as nextRows.
    ::Kokkos::View<PetscScalar*, MemorySpace> nextValues;

    /// @brief Matrix values.
    ::Kokkos::View<PetscScalar*, MemorySpace> values;

    /// @brief Coordinates sorted by row then column, as indices in the list.
    COOOrder<MemorySpace> order;

    /**
     * @brief Start in the sorted coordinates of each distinct coordinate,
     * and the list size.
     */
    ::Kokkos::View<Index*, MemorySpace> starts;

    /// @brief Sum of the values of each distinct coordinate.
    ::Kokkos::View<PetscScalar*, MemorySpace> merged;

    /// @brief Right-hand side.
    ::Kokkos::View<Real*, MemorySpace> rhs;

    /// @brief Indices of the linear integrator entries.
    ::Kokkos::View<Index*, MemorySpace> entryIndices;

    /// @brief Values of the linear integrator entries.
    ::Kokkos::View<Real*, MemorySpace> entryValues;

    /**
     * @brief Rows of the contributions to the right-hand side, listed in a
     * fixed order.
     */
    ::Kokkos::View<Index*, MemorySpace> contributionRows;

    /// @brief Contributions to the right-hand side.
    ::Kokkos::View<Real*, MemorySpace> contributionValues;

    /// @brief Number of listed contributions to the right-hand side.
    size_t contributions = 0;

    /// @brief Rows of the contributions grouped at the last assembly.
    ::Kokkos::View<Index*, MemorySpace> groupedRows;

    /**
     * @brief Positions of these contributions, grouped by row in list
     * order.
     */
    ::Kokkos::View<Index*, MemorySpace> groupedPositions;

    /**
     * @brief Start of every row group in groupedPositions, followed by the
     * contribution count.
     */
    ::Kokkos::View<Index*, MemorySpace> groupStarts;

    /// @brief Row of every group.
    ::Kokkos::View<Index*, MemorySpace> groupKeys;

    /**
     * @brief Linear integrators whose entries were grouped at the last
     * assembly, when the entries were the contributions themselves.
     */
    std::vector<COOSegment> vectorSegments;

    /// @brief Whether each degree of freedom is fixed.
    ::Kokkos::View<unsigned char*, MemorySpace> fixed;

    /// @brief Values of the fixed degrees of freedom.
    ::Kokkos::View<Real*, MemorySpace> fixedValues;

    /**
     * @brief Whether the right-hand side of each row is imposed, when
     * degrees of freedom are identified.
     */
    ::Kokkos::View<unsigned char*, MemorySpace> imposed;

    /// @brief Revision of the pattern, increased when the coordinates change.
    std::uint64_t pattern = 0;

    /// @brief Fixed degrees of freedom of the last assembly, on the host.
    std::vector<unsigned char> fixedMask;

    /// @brief Revision of fixedMask, increased when it changes.
    std::uint64_t fixedRevision = 0;

    /**
     * @brief Revisions of the pattern and of the fixed degrees of freedom
     * the eliminated entries were listed for.
     */
    std::pair<std::uint64_t, std::uint64_t> elimination{
      std::numeric_limits<std::uint64_t>::max(), 0};

    /// @brief Number of entries the eliminated entries were listed in.
    size_t eliminationCount = 0;

    /// @brief Entries of a free row and a fixed column, by increasing position.
    ::Kokkos::View<Index*, MemorySpace> movedEntries;

    /// @brief Row of each entry of movedEntries.
    ::Kokkos::View<Index*, MemorySpace> movedRows;

    /// @brief Column of each entry of movedEntries.
    ::Kokkos::View<Index*, MemorySpace> movedCols;

    /// @brief Entries of a fixed row.
    ::Kokkos::View<Index*, MemorySpace> clearedEntries;

    /// @brief Whether the matrix was preallocated with the expanded pattern.
    bool expanded = false;

    /// @brief Row coordinates expanded on the masters.
    ::Kokkos::View<Index*, MemorySpace> expandedRows;

    /// @brief Column coordinates expanded on the masters.
    ::Kokkos::View<Index*, MemorySpace> expandedCols;

    /// @brief Expanded coordinates being listed.
    ::Kokkos::View<Index*, MemorySpace> nextExpandedRows;

    /// @brief Expanded coordinates being listed.
    ::Kokkos::View<Index*, MemorySpace> nextExpandedCols;

    /// @brief Expanded matrix values.
    ::Kokkos::View<PetscScalar*, MemorySpace> expandedValues;

    /// @brief Expanded coordinates sorted by row then column.
    COOOrder<MemorySpace> expandedOrder;

    /// @brief Start of each distinct expanded coordinate, and the list size.
    ::Kokkos::View<Index*, MemorySpace> expandedStarts;

    /**
     * @brief Ghosted vector of two entries per row summing the fixed degrees
     * of freedom over the processes, kept while the distribution keeps its
     * layout.
     */
    std::shared_ptr<std::remove_pointer_t<::Vec>> exchange;

    /// @brief Distribution the exchange was created for.
    std::shared_ptr<const COODistribution> exchangeLayout;

    /**
     * @brief Local degrees of freedom held by several processes, the only
     * ones the exchange changes.
     */
    std::vector<Index> exchangeShared;

    /**
     * @brief Block of each degree of freedom of exchangeShared in the local
     * form of the exchange.
     */
    std::vector<PetscInt> exchangeBlocks;
  };

  /**
   * @brief PETSc linear system assembled in device memory, with value
   * constraints eliminated.
   *
   * The elimination is the one of SystemAssembly: rows of fixed degrees
   * of freedom become identity rows holding the fixed value, and the
   * columns of fixed degrees of freedom move to the right-hand side.
   * Eliminated entries keep their coordinates with a zero value, and every
   * diagonal entry is listed, so the pattern of the matrix does not depend
   * on the constraints.
   *
   * Integrators with the same layout between the same spaces share one
   * segment of the coordinate list, their values being summed in place, so
   * a sum of integrands such as the terms of an Oseen operator lists each
   * local matrix entry once.
   *
   * Identification constraints are resolved on the host, see
   * Constraints. The coordinate list is then expanded on the masters
   * of the rows and columns, the identified columns lifting their affine
   * values to the right-hand side, and the identified rows become the rows
   * of their constraints. The expanded pattern is kept as the plain one.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
}

#endif
