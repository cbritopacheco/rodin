/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_PETSC_H
#define RODIN_KOKKOS_ASSEMBLY_PETSC_H

/**
 * @file
 * @brief Device assembly of PETSc problems on resident Kokkos meshes.
 *
 * The integrators write their entries straight into coordinate lists in
 * device memory, which PETSc consumes through MatSetPreallocationCOO and
 * MatSetValuesCOO. The matrix and the vectors are created with the PETSc
 * Kokkos types, so a solve on the device never copies the system to the
 * host.
 */

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include "Rodin/Configure.h"
#include "Rodin/PETSc/Assembly/VectorSetup.h"
#ifdef RODIN_USE_MPI
#include "Rodin/PETSc/Assembly/ConstraintMap.h"
#endif
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/PETSc/Variational/BilinearForm.h"
#include "Rodin/PETSc/Variational/LinearForm.h"
#include "Rodin/PETSc/Variational/Problem.h"

#include "Rodin/Kokkos/Assembly/Kokkos.h"

#ifndef PETSC_HAVE_KOKKOS
#error "Device assembly of PETSc problems requires PETSc configured with Kokkos."
#endif

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Block of the coordinate list shared by integrators with the same
   * entries.
   */
  struct COOSegment
  {
    /// @brief Arrangement of the entries.
    ::Rodin::Kokkos::Variational::DeviceLayout layout;

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

    bool operator==(const COOSegment&) const = default;
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

    /// @brief Calls the function with the view that holds the indices.
    template <class Function>
    decltype(auto) visit(Function&& function) const
    {
      if (wide.extent(0) > 0)
        return function(wide);
      return function(narrow);
    }
  };

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
   * since the attributes of the cells may change. The buffers are reused
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

    COOCache(const COOCache&) {}

    COOCache& operator=(const COOCache& other)
    {
      if (this != &other)
        *this = COOCache();
      return *this;
    }

    COOCache(COOCache&&) = default;

    COOCache& operator=(COOCache&&) = default;

    /// @brief PETSc identifier of the matrix preallocated with the pattern, or -1.
    PetscInt64 matrix = -1;

    /// @brief PETSc identifier of the distributed right-hand side
    /// preallocated with the local coordinates, or -1.
    PetscInt64 vector = -1;

    /// @brief Whether the coordinates hold the pattern of the segments.
    bool indexed = false;

    /// @brief Segments of the pattern.
    std::vector<COOSegment> segments;

    /// @brief Row coordinates.
    ::Kokkos::View<Index*, MemorySpace> rows;

    /// @brief Column coordinates.
    ::Kokkos::View<Index*, MemorySpace> cols;

    /// @brief Row coordinates being listed.
    ::Kokkos::View<Index*, MemorySpace> nextRows;

    /// @brief Column coordinates being listed.
    ::Kokkos::View<Index*, MemorySpace> nextCols;

    /// @brief Matrix values.
    ::Kokkos::View<PetscScalar*, MemorySpace> values;

    /// @brief Coordinates sorted by row then column, as indices in the list.
    COOOrder<MemorySpace> order;

    /// @brief Start in the sorted coordinates of each distinct coordinate,
    /// and the list size.
    ::Kokkos::View<Index*, MemorySpace> starts;

    /// @brief Sum of the values of each distinct coordinate.
    ::Kokkos::View<PetscScalar*, MemorySpace> merged;

    /// @brief Right-hand side.
    ::Kokkos::View<Real*, MemorySpace> rhs;

    /// @brief Indices of the linear integrator entries.
    ::Kokkos::View<Index*, MemorySpace> entryIndices;

    /// @brief Values of the linear integrator entries.
    ::Kokkos::View<Real*, MemorySpace> entryValues;

    /// @brief Rows of the contributions to the right-hand side, listed in a
    /// fixed order.
    ::Kokkos::View<Index*, MemorySpace> contributionRows;

    /// @brief Contributions to the right-hand side.
    ::Kokkos::View<Real*, MemorySpace> contributionValues;

    /// @brief Number of listed contributions to the right-hand side.
    size_t contributions = 0;

    /// @brief Rows of the contributions grouped at the last assembly.
    ::Kokkos::View<Index*, MemorySpace> groupedRows;

    /// @brief Positions of these contributions, grouped by row in list
    /// order.
    ::Kokkos::View<Index*, MemorySpace> groupedPositions;

    /// @brief Start of every row group in groupedPositions, followed by the
    /// contribution count.
    ::Kokkos::View<Index*, MemorySpace> groupStarts;

    /// @brief Row of every group.
    ::Kokkos::View<Index*, MemorySpace> groupKeys;

    /// @brief Linear integrators whose entries were grouped at the last
    /// assembly, when the entries were the contributions themselves.
    std::vector<COOSegment> vectorSegments;

    /// @brief Whether each degree of freedom is fixed.
    ::Kokkos::View<unsigned char*, MemorySpace> fixed;

    /// @brief Values of the fixed degrees of freedom.
    ::Kokkos::View<Real*, MemorySpace> fixedValues;

    /// @brief Whether the right-hand side of each row is imposed, when
    /// degrees of freedom are identified.
    ::Kokkos::View<unsigned char*, MemorySpace> imposed;

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
  };

  /**
   * @brief PETSc linear system assembled in device memory, with value
   * constraints eliminated.
   *
   * The elimination is the one of DeviceSystem: rows of fixed degrees
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
   * DeviceConstraints. The coordinate list is then expanded on the masters
   * of the rows and columns, the identified columns lifting their affine
   * values to the right-hand side, and the identified rows become the rows
   * of their constraints. The expanded pattern is kept as the plain one.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  /**
   * @brief Distribution of the degrees of freedom of a system over MPI
   * processes.
   *
   * The local degrees of freedom, ghosts included, are those the device
   * integrators number. Each one has a global index and is owned by one
   * process, the one owning it in its space. The process holds the
   * contiguous range @f$ [b, e) @f$ of the rows of the PETSc layout, whose
   * length is its number of owned degrees of freedom. With a single field
   * this range is the owned range of the space. With several fields, each
   * numbered from its offset in the system as on the host, the owned
   * degrees of freedom of a field may lie in the rows of other processes,
   * and PETSc sends their entries there.
   */
  struct COODistribution
  {
    /// @brief Communicator of the processes.
    MPI_Comm comm = MPI_COMM_NULL;

    /// @brief First row of this process in the PETSc layout.
    PetscInt begin = 0;

    /// @brief One past the last row of this process in the PETSc layout.
    PetscInt end = 0;

    /// @brief Global number of degrees of freedom.
    PetscInt size = 0;

    /// @brief Global index of each local degree of freedom.
    std::vector<PetscInt> globals;

    /// @brief Whether this process owns each local degree of freedom.
    std::vector<unsigned char> owned;

    /// @brief Local index of the first degree of freedom of each field.
    std::vector<Index> starts;

    /// @brief Global index of the first degree of freedom of each field.
    std::vector<PetscInt> offsets;

    /// @brief Gets the number of owned degrees of freedom. @return Owned count.
    PetscInt getOwnedSize() const
    {
      return end - begin;
    }

    /**
     * @brief Appends the local degrees of freedom of a field.
     * @param fes Space numbering the local degrees of freedom of the field.
     * @param offset Global index of the first degree of freedom of the field.
     */
    template <class FES>
    void append(const FES& fes, PetscInt offset)
    {
      Index first = 0, last = 0;
      fes.getOwnershipRange(first, last);
      const size_t count = fes.getLocalSize();
      starts.push_back(static_cast<Index>(globals.size()));
      offsets.push_back(offset);
      globals.reserve(globals.size() + count);
      owned.reserve(owned.size() + count);
      for (size_t i = 0; i < count; ++i)
      {
        const Index g = fes.getGlobalIndex(static_cast<Index>(i));
        globals.push_back(offset + static_cast<PetscInt>(g));
        owned.push_back(g >= first && g < last);
      }
    }

    /**
     * @brief Gets the global offset of a field.
     * @param start Local index of the first degree of freedom of the field.
     * @return Global index of the first degree of freedom of the field.
     */
    PetscInt getOffset(Index start) const
    {
      const auto it = std::find(starts.begin(), starts.end(), start);
      assert(it != starts.end());
      return offsets[it - starts.begin()];
    }

    /// @brief Gets the local index of each global index held. @return Map from global to local indices.
    std::unordered_map<PetscInt, Index> getLocals() const
    {
      std::unordered_map<PetscInt, Index> res;
      res.reserve(globals.size());
      for (size_t i = 0; i < globals.size(); ++i)
        res.emplace(globals[i], static_cast<Index>(i));
      return res;
    }

    /**
     * @brief Lays out the rows of the processes in rank order; collective.
     * @param communicator Communicator of the processes.
     * @param total Global number of degrees of freedom.
     */
    void close(MPI_Comm communicator, PetscInt total)
    {
      comm = communicator;
      size = total;
      const PetscInt count =
        static_cast<PetscInt>(std::count(owned.begin(), owned.end(), 1));
      PetscInt first = 0;
      MPI_Exscan(&count, &first, 1, MPIU_INT, MPI_SUM, comm);
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      begin = rank == 0 ? 0 : first;
      end = begin + count;
    }

    /**
     * @brief Distribution of a distributed device space.
     * @param fes Space numbering the local degrees of freedom.
     * @return Distribution of the space.
     */
    template <class FES>
    static COODistribution of(const FES& fes)
    {
      COODistribution res;
      res.comm = fes.getMesh().getContext().getCommunicator();
      res.append(fes, 0);
      Index begin = 0, end = 0;
      fes.getOwnershipRange(begin, end);
      res.begin = static_cast<PetscInt>(begin);
      res.end = static_cast<PetscInt>(end);
      res.size = static_cast<PetscInt>(fes.getSize());
      return res;
    }
  };

  template <class Backend>
  class DeviceCOOSystem
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space of the kernels.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Memory space of the system.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Cache type.
      using CacheType = COOCache<MemorySpace>;

      /// @brief Constraints of the system.
      using ConstraintsType = DeviceConstraints<Backend>;

      /// @brief Expansions of the degrees of freedom on their masters.
      using Expansions = typename ConstraintsType::Expansions;

      /// @brief Device system whose grouping sums the right-hand side.
      using ReferenceType = DeviceSystem<Backend>;

      /// @brief Numbers of expanded contributions.
      using Sizes = typename ReferenceType::Sizes;

      /// @brief Contributions to the right-hand side reserved in the list.
      struct Contributions
      {
        /// @brief Rows of the listed contributions.
        ::Kokkos::View<Index*, MemorySpace> rows;

        /// @brief Listed contributions.
        ::Kokkos::View<Real*, MemorySpace> values;

        /// @brief Position of the first reserved contribution.
        size_t offset = 0;
      };

      /// @brief Device bilinear integrator.
      using BilinearType =
        ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>;

      /// @brief Device linear integrator.
      using LinearType =
        ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>;

      static_assert(std::is_same_v<PetscScalar, Real>,
        "Device assembly of PETSc problems requires real PETSc scalars.");

      /**
       * @brief Creates an empty system.
       * @param execution Execution space instance of the mesh.
       * @param rows Number of rows, the size of the test spaces.
       * @param cols Number of columns, the size of the trial spaces.
       */
      DeviceCOOSystem(const ExecutionSpace& execution, size_t rows, size_t cols)
        : m_execution(execution),
          m_rows(rows),
          m_cols(cols),
          m_constraints(rows, std::max(rows, cols))
      {
        if (std::max(rows, cols) >
          static_cast<size_t>(std::numeric_limits<PetscInt>::max()))
        {
          Alert::Exception() << "The system is too large for the PETSc index type."
                             << Alert::Raise;
        }
      }

      /**
       * @brief Assembles a Dirichlet condition and constrains its degrees of
       * freedom.
       * @param dbc Value or identification Dirichlet condition.
       * @param offset Offset of the constrained space in the system.
       * @param masterOffset Offset of the space of the masters of an
       * identification.
       */
      void fix(
        Variational::DirichletBCBase<Real>& dbc, size_t offset, size_t masterOffset)
      {
        m_constraints.fix(dbc, offset, masterOffset);
      }

      /**
       * @brief Assembles a Dirichlet condition whose masters, if any, are in
       * the constrained space.
       * @param dbc Value or identification Dirichlet condition.
       * @param offset Offset of the constrained space in the system.
       */
      void fix(Variational::DirichletBCBase<Real>& dbc, size_t offset)
      {
        fix(dbc, offset, offset);
      }

      /**
       * @brief Assembles a periodic condition and identifies its degrees of
       * freedom.
       * @param pbc Periodic condition.
       * @param offset Offset of the constrained space in the system.
       */
      void identify(Variational::PeriodicBCBase<Real>& pbc, size_t offset)
      {
        m_constraints.identify(pbc, offset);
      }

      /**
       * @brief Whether the matrix contributes to the right-hand side;
       * collective once the system is distributed, since a process may hold
       * degrees of freedom fixed by another one.
       */
      bool hasConstraints() const
      {
        const bool local = m_constraints.isConstrained() || m_constraints.isIdentified();
        return m_distribution ? isFixed(local) : local;
      }

      /**
       * @brief Distributes the system over MPI processes.
       * @param distribution Global numbering of the local rows and columns.
       *
       * The rows and columns of the system are then the local degrees of
       * freedom. Their entries are sent to the owners of their global
       * indices, which sum them. Every process eliminates the fixed degrees
       * of freedom it holds, and their owner alone writes the diagonal and
       * the fixed value.
       */
      void distribute(COODistribution distribution)
      {
        if (distribution.globals.size() != m_rows || (m_cols > 0 && m_rows != m_cols))
        {
          Alert::Exception()
            << "A distributed device system requires the local trial and test "
               "degrees of freedom to be the same."
            << Alert::Raise;
        }
        copy(m_owned, distribution.owned, "Rodin owned DOFs");
        m_distribution = std::move(distribution);
      }

      /**
       * @brief Distributes a system without constraints whose rows and
       * columns are numbered apart.
       * @param rows Global numbering of the local rows.
       * @param cols Global numbering of the local columns.
       */
      void distribute(COODistribution rows, COODistribution cols)
      {
        if (rows.globals.size() != m_rows || cols.globals.size() != m_cols)
        {
          Alert::Exception()
            << "The distributions do not number the local rows and columns."
            << Alert::Raise;
        }
        copy(m_owned, rows.owned, "Rodin owned DOFs");
        m_distribution = std::move(rows);
        m_columns = std::move(cols);
      }

      /// @brief Gets the distribution of a distributed system. @return Distribution.
      const COODistribution& getDistribution() const
      {
        assert(m_distribution);
        return *m_distribution;
      }

      /**
       * @brief Adds a local integrator to the system.
       * @param lbfi Local integrator, possibly negated or scaled by a real.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi,
        Index row, Index col)
      {
        const auto [integrator, sign] = getDeviceIntegrator<Backend>(lbfi);
        m_matrices.push_back({integrator, row, col, sign});
      }

      /**
       * @brief Adds a global integrator to the system.
       * @param gbfi Global integrator, possibly negated.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(const Variational::GlobalBilinearFormIntegratorBase<Real>& gbfi,
        Index row, Index col)
      {
        const auto [integrator, sign] = getDeviceIntegrator<Backend>(gbfi);
        m_matrices.push_back({integrator, row, col, sign});
      }

      /**
       * @brief Adds a linear integrator to the system.
       * @param lfi Integrator, possibly negated.
       * @param row Offset of its test space.
       * @param sign Factor of the contributions.
       */
      void addVector(
        const Variational::LinearFormIntegratorBase<Real>& lfi, Index row, Real sign)
      {
        const auto [integrator, s] = getDeviceIntegrator<Backend>(lfi, sign);
        m_vectors.push_back({integrator, row, s});
      }

      /**
       * @brief Adds the entries of a preassembled operator to the matrix.
       * @param op Operator, alive until finish() returns.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(::Mat op, Index row, Index col)
      {
        m_operators.push_back({op, row, col});
      }

      /**
       * @brief Adds a preassembled vector to the right-hand side.
       * @param vec Vector, alive until finish() returns.
       * @param row Offset of its test space.
       */
      void addVector(::Vec vec, Index row)
      {
        m_preassembled.push_back({vec, row});
      }

      /**
       * @brief Integrates, eliminates the constraints and fills the PETSc
       * objects.
       * @param cache Buffers and pattern of the previous assemblies.
       * @param A Matrix, filled when not null.
       * @param b Right-hand side, filled when not null.
       */
      void finish(CacheType& cache, ::Mat A, ::Vec b)
      {
        typename ConstraintsType::Resolved constraints;
        if (m_distribution && isFixed(m_constraints.isIdentified()))
          constraints = resolveDistributed();
        else if (m_constraints.isIdentified())
          constraints = m_constraints.resolve(m_execution);
        const size_t n = std::max(m_rows, m_cols);
        const bool expanded = constraints.isExpanded();
        // A degree of freedom fixed on one process is fixed on all of them.
        const bool fixed = m_distribution ? isFixed(m_constraints.isConstrained())
                                          : m_constraints.isConstrained();
        if (m_distribution && fixed)
        {
          auto mask = m_constraints.getFixed();
          auto values = m_constraints.getValues();
          mask.resize(m_rows, 0);
          values.resize(m_rows, Real(0));
          synchronize(mask, values);
          copy(cache.fixed, mask, "Rodin fixed DOFs");
          copy(cache.fixedValues, values, "Rodin fixed values");
        }
        else if (fixed || expanded)
        {
          copy(cache.fixed, m_constraints.getFixed(), "Rodin fixed DOFs");
          copy(cache.fixedValues, m_constraints.getValues(), "Rodin fixed values");
        }
        if (expanded)
          copy(cache.imposed, constraints.imposed, "Rodin imposed DOFs");
        resize(cache.rhs, n, "Rodin right-hand side");
        ::Kokkos::deep_copy(m_execution, cache.rhs, Real(0));
        cache.contributions = 0;

        if (A || !m_matrices.empty() || !m_operators.empty())
        {
          // The constraints are eliminated after the expansion.
          integrateMatrices(cache, fixed && !expanded);
          if (expanded)
            expandMatrices(cache, constraints, fixed);
        }
        integrateVectors(cache, fixed, expanded ? &constraints.expansions : nullptr);
        sumContributions(cache);
        if (expanded)
          imposeFixed(
            m_execution, m_rows, cache.imposed, cache.fixedValues, m_owned, cache.rhs);
        else if (fixed)
        {
          imposeFixed(
            m_execution, m_rows, cache.fixed, cache.fixedValues, m_owned, cache.rhs);
        }
        m_execution.fence("Rodin PETSc system");

        if (A && expanded)
        {
          setMatrix(cache, A, true, cache.expandedRows, cache.expandedCols,
            cache.expandedOrder, cache.expandedStarts, cache.expandedValues);
        }
        else if (A)
        {
          setMatrix(cache, A, false, cache.rows, cache.cols, cache.order, cache.starts,
            cache.values);
        }
        if (b && m_distribution)
          setDistributedVector(cache, b);
        else if (b)
          setVector(b, cache.rhs);
        // PETSc reads the buffers on its own execution space instance, and
        // the next assembly overwrites them.
        ::Kokkos::fence("Rodin PETSc system");
      }

      /// @brief Adds the block offsets to the coordinates of entry e.
      template <class Indices>
      struct Shift
      {
        Indices rows;
        Indices cols;
        Index row;
        Index col;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          rows(e) += row;
          cols(e) += col;
        }
      };

      /// @brief Adds the block offsets to the coordinates of a segment.
      template <class Indices>
      static void shift(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, size_t offset, size_t count, Index row, Index col)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::shift",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, offset, offset + count),
          Shift<Indices>{ rows, cols, row, col });
      }

      /// @brief Lists the diagonal entry i after the given offset.
      template <class Indices>
      struct DiagonalList
      {
        size_t offset;
        Indices rows;
        Indices cols;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          rows(offset + i) = i;
          cols(offset + i) = i;
        }
      };

      /// @brief Lists every diagonal entry.
      template <class Indices>
      static void listDiagonal(const ExecutionSpace& execution, size_t offset,
        size_t count, const Indices& rows, const Indices& cols)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::diagonal",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          DiagonalList<Indices>{ offset, rows, cols });
      }

      /// @brief Writes the diagonal value of row i, one when it is fixed.
      template <class Fixed, class Owned, class Scalars>
      struct DiagonalFill
      {
        size_t offset;
        bool hasFixed;
        Fixed fixed;
        Owned owned;
        Scalars coo;
        bool distributed;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          const bool one = hasFixed && fixed(i) && (!distributed || owned(i));
          coo(offset + i) = one ? PetscScalar(1) : PetscScalar(0);
        }
      };

      /// @brief Writes the diagonal values, one on fixed rows.
      template <class Fixed, class Owned, class Scalars>
      static void fillDiagonal(const ExecutionSpace& execution, size_t offset,
        size_t count, bool hasFixed, const Fixed& fixed, const Owned& owned,
        const Scalars& coo)
      {
        const bool distributed = owned.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::diagonal",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          DiagonalFill<Fixed, Owned, Scalars>{
            offset, hasFixed, fixed, owned, coo, distributed });
      }

      /// @brief Counts entry e when its value moves to the right-hand side.
      template <class Indices, class Fixed>
      struct EliminatedCount
      {
        Indices rows;
        Indices cols;
        Fixed fixed;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += !fixed(rows(e)) && fixed(cols(e));
        }
      };

      /**
       * @brief Counts the matrix values moved to the right-hand side by the
       * columns of the fixed degrees of freedom.
       */
      template <class Indices, class Fixed>
      static size_t countEliminated(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Fixed& fixed)
      {
        size_t res = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::countEliminated",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          EliminatedCount<Indices, Fixed>{ rows, cols, fixed }, res);
        return res;
      }

      /**
       * @brief Clears entry e when its row is fixed and lists its value when
       * its column alone is fixed, at the scanned position.
       */
      template <class Indices, class Fixed, class Values, class Scalars>
      struct Elimination
      {
        Indices rows;
        Indices cols;
        Fixed fixed;
        Values values;
        Scalars coo;
        ::Kokkos::View<Index*, MemorySpace> resRows;
        ::Kokkos::View<Real*, MemorySpace> resValues;
        size_t offset;

        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, size_t& update, const bool final) const
        {
          const Index r = rows(e);
          const Index c = cols(e);
          const bool moved = !fixed(r) && fixed(c);
          if (final && (moved || fixed(r)))
          {
            if (moved)
            {
              resRows(offset + update) = r;
              resValues(offset + update) = -(coo(e) * values(c));
            }
            coo(e) = 0;
          }
          update += moved;
        }
      };

      /**
       * @brief Clears the rows of the fixed degrees of freedom and lists the
       * values of their columns moved to the right-hand side.
       * @param res Receives the moved values, in coordinate order.
       */
      template <class Indices, class Fixed, class Values, class Scalars>
      static void eliminate(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Fixed& fixed,
        const Values& values, const Scalars& coo, const Contributions& res)
      {
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::PETSc::eliminate",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          Elimination<Indices, Fixed, Values, Scalars>{
            rows, cols, fixed, values, coo, res.rows, res.values, res.offset });
      }

      /**
       * @brief Number of contributions of an entry of a linear integrator.
       * @param row Offset of the space of the integrator.
       * @param x Expansions on the masters, used when @p expanded is set.
       * @param expanded Whether an entry goes to the masters of its row, with
       * their coefficients, rather than to its row when it is free.
       */
      template <class Indices, class Values, class Fixed>
      KOKKOS_INLINE_FUNCTION static size_t countEntry(size_t e, const Indices& indices,
        const Values& values, Index row, bool hasFixed, const Fixed& fixed,
        const Expansions& x, bool expanded)
      {
        if (values(e) == Real(0))
          return 0;
        const Index i = row + indices(e);
        if (expanded)
          return x.count(i);
        return !hasFixed || !fixed(i);
      }

      /// @brief Counts the contributions of entry e of a linear integrator.
      template <class Indices, class Values, class Fixed>
      struct VectorCount
      {
        Indices indices;
        Values values;
        Index row;
        bool hasFixed;
        Fixed fixed;
        Expansions x;
        bool expanded;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += countEntry(e, indices, values, row, hasFixed, fixed, x, expanded);
        }
      };

      /**
       * @brief Counts the contributions of the entries of a linear
       * integrator to the right-hand side.
       * @param begin Position of the first entry. @param count Entry count.
       *
       * See countEntry for the other parameters.
       */
      template <class Indices, class Values, class Fixed>
      static size_t countVector(const ExecutionSpace& execution, size_t begin,
        size_t count, const Indices& indices, const Values& values, Index row,
        bool hasFixed, const Fixed& fixed, const Expansions& x, bool expanded)
      {
        size_t res = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::countVector",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, begin, begin + count),
          VectorCount<Indices, Values, Fixed>{
            indices, values, row, hasFixed, fixed, x, expanded },
          res);
        return res;
      }

      /**
       * @brief Lists the contributions of entry e of a linear integrator at
       * the scanned position.
       */
      template <class Indices, class Values, class Fixed>
      struct VectorList
      {
        Indices indices;
        Values values;
        Index row;
        bool hasFixed;
        Fixed fixed;
        Expansions x;
        bool expanded;
        ::Kokkos::View<Index*, MemorySpace> resRows;
        ::Kokkos::View<Real*, MemorySpace> resValues;
        size_t offset;

        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, size_t& update, const bool final) const
        {
          const size_t n =
            countEntry(e, indices, values, row, hasFixed, fixed, x, expanded);
          if (final && n > 0)
          {
            const Real v = values(e);
            const Index i = row + indices(e);
            if (expanded)
            {
              const Index r0 = x.starts(i);
              for (size_t r = 0; r < n; ++r)
              {
                resRows(offset + update + r) = x.indices(r0 + r);
                resValues(offset + update + r) = x.coefficients(r0 + r) * v;
              }
            }
            else
            {
              resRows(offset + update) = i;
              resValues(offset + update) = v;
            }
          }
          update += n;
        }
      };

      /**
       * @brief Lists the nonzero entries of a linear integrator as
       * contributions to the right-hand side.
       * @param res Receives the contributions, counted by countVector.
       *
       * The contributions follow the order of the entries, and of the masters
       * of an entry. See countVector for the other parameters.
       */
      template <class Indices, class Values, class Fixed>
      static void listVector(const ExecutionSpace& execution, size_t begin, size_t count,
        const Indices& indices, const Values& values, Index row, bool hasFixed,
        const Fixed& fixed, const Expansions& x, bool expanded, const Contributions& res)
      {
        const auto policy =
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, begin, begin + count);
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::PETSc::vector", policy,
          VectorList<Indices, Values, Fixed>{ indices, values, row, hasFixed, fixed, x,
            expanded, res.rows, res.values, res.offset });
      }

      /// @brief Writes the fixed value of row i when it is fixed.
      template <class Fixed, class Values, class Owned, class RHS>
      struct FixedImposition
      {
        Fixed fixed;
        Values values;
        Owned owned;
        RHS rhs;
        bool distributed;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          if (fixed(i))
            rhs(i) = !distributed || owned(i) ? values(i) : Real(0);
        }
      };

      /// @brief Writes the fixed values on the fixed rows.
      template <class Fixed, class Values, class Owned, class RHS>
      static void imposeFixed(const ExecutionSpace& execution, size_t count,
        const Fixed& fixed, const Values& values, const Owned& owned, const RHS& rhs)
      {
        const bool distributed = owned.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::fixed",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          FixedImposition<Fixed, Values, Owned, RHS>{
            fixed, values, owned, rhs, distributed });
      }

      /// @brief Writes the row of entry e of a linear integrator.
      template <class Indices>
      struct EntryRows
      {
        Indices indices;
        Index row;
        Indices res;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          res(e) = row + indices(e);
        }
      };

      /**
       * @brief Writes the rows of the entries of a linear integrator.
       * @param indices Entry indices. @param begin Position of the first
       * entry. @param end Position past the last one.
       * @param row Offset of the test space. @param res Receives the rows.
       */
      template <class Indices>
      static void shiftEntries(const ExecutionSpace& execution, const Indices& indices,
        size_t begin, size_t end, Index row, const Indices& res)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::rows",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, begin, end),
          EntryRows<Indices>{ indices, row, res });
      }

      /// @brief Numbers entry i of a dense vector.
      template <class Indices>
      struct Sequence
      {
        Indices res;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const { res(i) = i; }
      };

      /// @brief Numbers the entries of a dense vector.
      template <class Indices>
      static void sequence(const ExecutionSpace& execution, const Indices& res)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::sequence",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, res.extent(0)),
          Sequence<Indices>{ res });
      }

      /**
       * @brief Counts the contributions of entry e expanded on the masters,
       * and the affine values it lifts.
       */
      template <class Indices, class Scalars>
      struct ExpandedCount
      {
        Indices rows;
        Indices cols;
        Scalars values;
        Expansions x;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, Sizes& sum) const
        {
          if (values(e) == PetscScalar(0))
            return;
          const size_t nr = x.count(rows(e));
          Sizes c;
          c.entries = nr * x.count(cols(e));
          c.lifts = x.values(cols(e)) != Real(0) ? nr : 0;
          sum += c;
        }
      };

      /**
       * @brief Counts the contributions of the coordinate list expanded on
       * the masters, and the affine values they lift.
       */
      template <class Indices, class Scalars>
      static Sizes countExpanded(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Scalars& values,
        const Expansions& x)
      {
        Sizes res;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::countExpanded",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          ExpandedCount<Indices, Scalars>{ rows, cols, values, x }, res);
        return res;
      }

      /**
       * @brief Expands entry e on the masters and lists the affine values it
       * lifts, at the scanned positions.
       */
      template <class Indices, class Scalars>
      struct Expansion
      {
        Indices rows;
        Indices cols;
        Scalars values;
        Expansions x;
        Indices resRows;
        Indices resCols;
        Scalars resValues;
        ::Kokkos::View<Index*, MemorySpace> liftRows;
        ::Kokkos::View<Real*, MemorySpace> liftValues;
        size_t offset;

        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, Sizes& update, const bool final) const
        {
          const PetscScalar v = values(e);
          if (v == PetscScalar(0))
            return;
          const size_t nr = x.count(rows(e));
          const size_t nc = x.count(cols(e));
          const Real d = x.values(cols(e));
          if (final)
          {
            const Index r0 = x.starts(rows(e));
            const Index c0 = x.starts(cols(e));
            for (size_t ri = 0; ri < nr; ++ri)
            {
              const Real rv = x.coefficients(r0 + ri) * v;
              for (size_t ci = 0; ci < nc; ++ci)
              {
                const size_t q = update.entries + ri * nc + ci;
                resRows(q) = x.indices(r0 + ri);
                resCols(q) = x.indices(c0 + ci);
                resValues(q) = rv * x.coefficients(c0 + ci);
              }
              if (d != Real(0))
              {
                liftRows(offset + update.lifts + ri) = x.indices(r0 + ri);
                liftValues(offset + update.lifts + ri) = -(rv * d);
              }
            }
          }
          Sizes c;
          c.entries = nr * nc;
          c.lifts = d != Real(0) ? nr : 0;
          update += c;
        }
      };

      /**
       * @brief Expands the coordinate list on the masters.
       * @param lifts Receives the lifted affine values, in coordinate order.
       *
       * The nonzero value @f$ a @f$ of row @f$ i @f$ and column @f$ j @f$
       * becomes the values @f$ (c_r a) c_s @f$ of the masters @f$ r @f$ of
       * @f$ i @f$ and @f$ s @f$ of @f$ j @f$, and lifts @f$ -(c_r a) d_j @f$
       * to row @f$ r @f$ when @f$ j @f$ has a nonzero affine value
       * @f$ d_j @f$, as the host assemblies do.
       */
      template <class Indices, class Scalars>
      static void expand(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Scalars& values,
        const Expansions& x, const Indices& resRows, const Indices& resCols,
        const Scalars& resValues, const Contributions& lifts)
      {
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::PETSc::expand",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          Expansion<Indices, Scalars>{ rows, cols, values, x, resRows, resCols,
            resValues, lifts.rows, lifts.values, lifts.offset });
      }

      /// @brief Counts position i when two index lists differ there.
      template <class Indices, class Previous>
      struct DifferenceCount
      {
        Indices indices;
        Previous previous;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t i, size_t& sum) const
        {
          sum += indices(i) != previous(i);
        }
      };

      /// @brief Counts the positions at which two index lists differ.
      template <class Indices, class Previous>
      static size_t countDifferences(
        const ExecutionSpace& execution, const Indices& indices, const Previous& previous)
      {
        size_t differences = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::differences",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, indices.extent(0)),
          DifferenceCount<Indices, Previous>{ indices, previous }, differences);
        return differences;
      }

      /// @brief Counts coordinate e when it differs between two patterns.
      template <class Indices>
      struct ChangeCount
      {
        Indices rows;
        Indices cols;
        Indices previousRows;
        Indices previousCols;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += rows(e) != previousRows(e) || cols(e) != previousCols(e);
        }
      };

      /// @brief Counts the coordinates that differ between two patterns.
      template <class Indices>
      static size_t countChanges(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, const Indices& previousRows, const Indices& previousCols)
      {
        size_t changes = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::pattern",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, rows.extent(0)),
          ChangeCount<Indices>{ rows, cols, previousRows, previousCols }, changes);
        return changes;
      }

      /**
       * @brief Sorts the coordinates by row then column.
       * @param columns Column count of the system.
       * @param segments Segments of the coordinates.
       * @param offsets Position of the first entry of every segment, the
       * entries past the last segment being stored in rank order.
       * @param rows Row coordinates. @param cols Column coordinates.
       * @param order Receives the indices of the sorted coordinates.
       *
       * The coordinates are handed to the sort in rank order, so that the
       * duplicates of a coordinate follow one another as the integrators
       * ranked them, whatever the spacing of their entries.
       */
      template <class Indices>
      static void sortPattern(const ExecutionSpace& execution, size_t columns,
        const std::vector<COOSegment>& segments, const std::vector<size_t>& offsets,
        const Indices& rows, const Indices& cols, COOOrder<MemorySpace>& order)
      {
        const size_t count = rows.extent(0);
        order = {};
        if (count <= std::numeric_limits<std::uint32_t>::max())
        {
          order.narrow = decltype(order.narrow)(
            ::Kokkos::view_alloc(
              execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO order")),
            count);
          sortPattern(execution, columns, segments, offsets, rows, cols, order.narrow);
        }
        else
        {
          order.wide = decltype(order.wide)(
            ::Kokkos::view_alloc(
              execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO order")),
            count);
          sortPattern(execution, columns, segments, offsets, rows, cols, order.wide);
        }
      }

      /**
       * @brief Writes the sort key and the index of the coordinate of rank r
       * in a segment.
       */
      template <class Indices, class Order>
      struct PatternKeys
      {
        using Position = typename Order::value_type;

        size_t begin;
        size_t size;
        size_t spacing;
        ::Kokkos::View<std::uint64_t*, MemorySpace> keys;
        Indices rows;
        Indices cols;
        size_t columns;
        Order order;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t r) const
        {
          const size_t e =
            begin + ::Rodin::Kokkos::Variational::getEntryPosition(r, size, spacing);
          keys(begin + r) = std::uint64_t(rows(e)) * columns + cols(e);
          order(begin + r) = static_cast<Position>(e);
        }
      };

      /// @brief Sorts the coordinates into the given indices.
      template <class Indices, class Order>
      static void sortPattern(const ExecutionSpace& execution, size_t columns,
        const std::vector<COOSegment>& segments, const std::vector<size_t>& offsets,
        const Indices& rows, const Indices& cols, const Order& order)
      {
        const size_t count = rows.extent(0);
        if (count == 0)
          return;
        const ::Kokkos::View<std::uint64_t*, MemorySpace> keys(
          ::Kokkos::view_alloc(
            execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO keys")),
          count);
        for (size_t s = 0; s <= segments.size(); ++s)
        {
          const size_t begin = offsets[s];
          const size_t size = (s < segments.size() ? offsets[s + 1] : count) - begin;
          const size_t spacing = s < segments.size() ? segments[s].spacing : 1;
          ::Kokkos::parallel_for(
            "Rodin::Assembly::Kokkos::PETSc::keys",
            ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, size),
            PatternKeys<Indices, Order>{
              begin, size, spacing, keys, rows, cols, columns, order });
        }
        ::Kokkos::Experimental::sort_by_key(execution, keys, order);
      }

      /**
       * @brief Whether the sorted coordinate i differs from the previous one.
       */
      template <class Indices, class Order>
      struct DistinctHead
      {
        Indices rows;
        Indices cols;
        Order order;

        KOKKOS_INLINE_FUNCTION bool operator()(const size_t i) const
        {
          const Index e = order(i);
          if (i == 0)
            return true;
          const Index p = order(i - 1);
          return rows(e) != rows(p) || cols(e) != cols(p);
        }
      };

      /// @brief Counts the sorted coordinate i when it starts a distinct one.
      template <class Indices, class Order>
      struct DistinctCount
      {
        DistinctHead<Indices, Order> head;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t i, size_t& sum) const
        {
          sum += head(i);
        }
      };

      /**
       * @brief Writes the sorted coordinate i at the scanned position when it
       * starts a distinct one.
       */
      template <class Indices, class Order>
      struct DistinctStartList
      {
        DistinctHead<Indices, Order> head;
        Indices starts;

        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t i, size_t& update, const bool final) const
        {
          if (head(i))
          {
            if (final)
              starts(update) = i;
            ++update;
          }
        }
      };

      /**
       * @brief Finds where each distinct coordinate starts in the sorted
       * coordinates.
       * @returns Starts of the distinct coordinates, followed by the list size.
       */
      template <class Indices, class Order>
      static Indices distinctStarts(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, const Order& order)
      {
        const size_t count = order.extent(0);
        const DistinctHead<Indices, Order> head{ rows, cols, order };
        size_t distinct = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::distinct",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          DistinctCount<Indices, Order>{ head }, distinct);
        Indices starts(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                         std::string("Rodin COO starts")),
          distinct + 1);
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::PETSc::starts",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          DistinctStartList<Indices, Order>{ head, starts });
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(starts, distinct), count);
        return starts;
      }

      /// @brief Lists the distinct coordinate u in the PETSc index type.
      template <class Indices, class Order, class Coordinates>
      struct DistinctList
      {
        Indices rows;
        Indices cols;
        Order order;
        Indices starts;
        Coordinates distinctRows;
        Coordinates distinctCols;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          const Index e = order(starts(u));
          distinctRows(u) = static_cast<PetscInt>(rows(e));
          distinctCols(u) = static_cast<PetscInt>(cols(e));
        }
      };

      /// @brief Lists the distinct coordinates in the PETSc index type.
      template <class Indices, class Order, class Coordinates>
      static void listDistinct(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, const Order& order, const Indices& starts,
        const Coordinates& distinctRows, const Coordinates& distinctCols)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::distinctCoordinates",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, distinctRows.extent(0)),
          DistinctList<Indices, Order, Coordinates>{
            rows, cols, order, starts, distinctRows, distinctCols });
      }

      /// @brief Sums the values of the distinct coordinate u.
      template <class Order, class Indices, class Scalars>
      struct Merge
      {
        Order order;
        Indices starts;
        Scalars values;
        Scalars merged;

        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          PetscScalar sum = 0;
          for (Index i = starts(u); i < starts(u + 1); ++i)
            sum += values(order(i));
          merged(u) = sum;
        }
      };

      /// @brief Sums the values of each distinct coordinate.
      template <class Order, class Indices, class Scalars>
      static void merge(const ExecutionSpace& execution, const Order& order,
        const Indices& starts, const Scalars& values, const Scalars& merged)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::merge",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, merged.extent(0)),
          Merge<Order, Indices, Scalars>{ order, starts, values, merged });
      }

    private:
      struct Matrix
      {
        const BilinearType* integrator;
        Index row;
        Index col;
        Real sign;
      };

      struct Vector
      {
        const LinearType* integrator;
        Index row;
        Real sign;
      };

      struct Operator
      {
        ::Mat op;
        Index row;
        Index col;
      };

      struct Preassembled
      {
        ::Vec vec;
        Index row;
      };

      /// @brief Coordinates and values of the preassembled operators, by row.
      struct Entries
      {
        std::vector<Index> rows;
        std::vector<Index> cols;
        std::vector<PetscScalar> values;
      };

      /**
       * @brief Lists the entries of the preassembled operators.
       *
       * A distributed operator lists its owned rows, which are complete, in
       * the global numbering of its spaces. Their entries are numbered
       * locally, as the integrated ones, so that they are eliminated and sent
       * alike.
       */
      Entries listOperators() const
      {
        Entries res;
        std::unordered_map<PetscInt, Index> locals;
        if (m_distribution && !m_operators.empty())
          locals = m_distribution->getLocals();
        const auto local = [&](PetscInt global) {
          const auto it = locals.find(global);
          if (it == locals.end())
          {
            Alert::Exception()
              << "A preassembled operator couples a degree of freedom the "
                 "process does not hold."
              << Alert::Raise;
          }
          return it->second;
        };
        for (const auto& o : m_operators)
        {
          PetscInt begin = 0, end = 0;
          PetscErrorCode ierr = MatGetOwnershipRange(o.op, &begin, &end);
          assert(ierr == PETSC_SUCCESS);
          const PetscInt rowOffset = m_distribution ? m_distribution->getOffset(o.row) : 0;
          const PetscInt colOffset = m_distribution ? m_distribution->getOffset(o.col) : 0;
          for (PetscInt i = begin; i < end; ++i)
          {
            PetscInt count = 0;
            const PetscInt* cols = nullptr;
            const PetscScalar* values = nullptr;
            ierr = MatGetRow(o.op, i, &count, &cols, &values);
            assert(ierr == PETSC_SUCCESS);
            for (PetscInt k = 0; k < count; ++k)
            {
              if (m_distribution)
              {
                res.rows.push_back(local(rowOffset + i));
                res.cols.push_back(local(colOffset + cols[k]));
              }
              else
              {
                res.rows.push_back(o.row + static_cast<Index>(i));
                res.cols.push_back(o.col + static_cast<Index>(cols[k]));
              }
              res.values.push_back(values[k]);
            }
            ierr = MatRestoreRow(o.op, i, &count, &cols, &values);
            assert(ierr == PETSC_SUCCESS);
          }
          (void)ierr;
        }
        return res;
      }

      /// @brief Unmanaged host view of a vector.
      template <class T>
      static auto host(const std::vector<T>& data)
      {
        return ::Kokkos::View<const T*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(data.data(), data.size());
      }

      template <class T>
      void resize(
        ::Kokkos::View<T*, MemorySpace>& view, size_t size, const char* label) const
      {
        if (view.extent(0) != size)
        {
          view = ::Kokkos::View<T*, MemorySpace>(
            ::Kokkos::view_alloc(
              m_execution, ::Kokkos::WithoutInitializing, std::string(label)),
            size);
        }
      }

      template <class T>
      void copy(::Kokkos::View<T*, MemorySpace>& view, const std::vector<T>& host,
        const char* label) const
      {
        resize(view, host.size(), label);
        const ::Kokkos::View<const T*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>
          source(host.data(), host.size());
        ::Kokkos::deep_copy(m_execution, view, source);
      }

      /**
       * @brief Integrates the local integrators into the coordinate list.
       */
      void integrateMatrices(CacheType& cache, bool fixed)
      {
        // Integrators with equal layouts and blocks share their entries,
        // listed on the diagonals of the component blocks when every one of
        // them vanishes elsewhere.
        std::vector<COOSegment> blocks;
        std::vector<bool> diagonals;
        std::vector<size_t> block(m_matrices.size());
        for (size_t k = 0; k < m_matrices.size(); ++k)
        {
          const auto& m = m_matrices[k];
          COOSegment key{m.integrator->getLayout(), m.row, m.col};
          const auto it = std::find(blocks.begin(), blocks.end(), key);
          block[k] = it - blocks.begin();
          if (it == blocks.end())
          {
            blocks.push_back(std::move(key));
            diagonals.push_back(true);
          }
          diagonals[block[k]] = diagonals[block[k]] && m.integrator->isBlockDiagonal();
        }
        std::vector<COOSegment> segments;
        std::vector<size_t> segment(m_matrices.size());
        bool restricted = false;
        for (size_t k = 0; k < m_matrices.size(); ++k)
        {
          const auto& m = m_matrices[k];
          const bool diagonal = diagonals[block[k]];
          const size_t entries = m.integrator->prepare(diagonal);
          COOSegment key{m.integrator->getLayout(), m.row, m.col, entries,
            m.integrator->getSpacing(), diagonal};
          restricted = restricted || !key.layout.attributes.empty() ||
            !key.layout.trialAttributes.empty();
          const auto it = std::find(segments.begin(), segments.end(), key);
          segment[k] = it - segments.begin();
          if (it == segments.end())
            segments.push_back(std::move(key));
        }
        std::vector<size_t> offsets(segments.size() + 1, 0);
        for (size_t s = 0; s < segments.size(); ++s)
          offsets[s + 1] = offsets[s] + segments[s].count;
        // The preassembled entries precede the diagonal, after the segments.
        const Entries entries = listOperators();
        const size_t tail = offsets.back() + entries.values.size();
        const size_t diagonal = std::min(m_rows, m_cols);
        const size_t count = tail + diagonal;

        const bool index = !cache.indexed || restricted || !m_operators.empty() ||
          cache.segments != segments || cache.rows.extent(0) != count;
        resize(cache.values, count, "Rodin COO values");
        if (index)
        {
          resize(cache.nextRows, count, "Rodin COO rows");
          resize(cache.nextCols, count, "Rodin COO cols");
        }
        const ::Kokkos::View<Index*, MemorySpace> none;
        std::vector<bool> filled(segments.size(), false);
        ::Rodin::Kokkos::Internal::clearLaunches();
        for (size_t k = 0; k < m_matrices.size(); ++k)
        {
          const auto& m = m_matrices[k];
          const size_t s = segment[k];
          const bool first = !filled[s];
          m.integrator->fill(cache.values, offsets[s], m.sign, !first,
            first && index ? cache.nextRows : none,
            first && index ? cache.nextCols : none);
          ::Rodin::Kokkos::Internal::checkLaunches("A device bilinear integrator");
          filled[s] = true;
        }
        const auto range = std::make_pair(offsets.back(), tail);
        ::Kokkos::deep_copy(
          m_execution, ::Kokkos::subview(cache.values, range), host(entries.values));
        if (index)
        {
          for (size_t s = 0; s < segments.size(); ++s)
          {
            shift(m_execution, cache.nextRows, cache.nextCols, offsets[s],
              segments[s].count, segments[s].row, segments[s].col);
          }
          ::Kokkos::deep_copy(
            m_execution, ::Kokkos::subview(cache.nextRows, range), host(entries.rows));
          ::Kokkos::deep_copy(
            m_execution, ::Kokkos::subview(cache.nextCols, range), host(entries.cols));
          listDiagonal(m_execution, tail, diagonal, cache.nextRows, cache.nextCols);
          const bool changed = !cache.indexed || cache.rows.extent(0) != count ||
            countChanges(
              m_execution, cache.nextRows, cache.nextCols, cache.rows, cache.cols) > 0;
          std::swap(cache.rows, cache.nextRows);
          std::swap(cache.cols, cache.nextCols);
          cache.segments = std::move(segments);
          cache.indexed = true;
          if (changed)
          {
            sortPattern(m_execution, m_cols, cache.segments, offsets, cache.rows,
              cache.cols, cache.order);
            cache.starts = cache.order.visit([&](const auto& order) {
              return distinctStarts(m_execution, cache.rows, cache.cols, order);
            });
            cache.matrix = -1;
          }
        }
        fillDiagonal(
          m_execution, tail, diagonal, fixed, cache.fixed, m_owned, cache.values);
        if (fixed)
        {
          const size_t moved =
            countEliminated(m_execution, tail, cache.rows, cache.cols, cache.fixed);
          eliminate(m_execution, tail, cache.rows, cache.cols, cache.fixed,
            cache.fixedValues, cache.values, reserve(cache, moved));
        }
      }

      /**
       * @brief Expands the coordinate list on the masters, appends the rows
       * of the identified degrees of freedom and every diagonal entry, then
       * eliminates the fixed degrees of freedom.
       */
      void expandMatrices(CacheType& cache,
        const typename ConstraintsType::Resolved& constraints, bool fixed)
      {
        const size_t listed = cache.values.extent(0);
        const Sizes sizes = countExpanded(m_execution, listed, cache.rows, cache.cols,
          cache.values, constraints.expansions);
        const size_t entries = sizes.entries;
        const size_t extra = constraints.rows.size();
        const size_t diagonal = std::min(m_rows, m_cols);
        const size_t count = entries + extra + diagonal;
        resize(cache.nextExpandedRows, count, "Rodin expanded COO rows");
        resize(cache.nextExpandedCols, count, "Rodin expanded COO cols");
        resize(cache.expandedValues, count, "Rodin expanded COO values");
        expand(m_execution, listed, cache.rows, cache.cols, cache.values,
          constraints.expansions, cache.nextExpandedRows, cache.nextExpandedCols,
          cache.expandedValues, reserve(cache, sizes.lifts));
        const auto range = std::make_pair(entries, entries + extra);
        ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(cache.nextExpandedRows, range),
          host(constraints.rows));
        ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(cache.nextExpandedCols, range),
          host(constraints.cols));
        ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(cache.expandedValues, range),
          host(constraints.values));
        listDiagonal(m_execution, entries + extra, diagonal, cache.nextExpandedRows,
          cache.nextExpandedCols);
        const bool changed = cache.expandedRows.extent(0) != count ||
          countChanges(m_execution, cache.nextExpandedRows, cache.nextExpandedCols,
            cache.expandedRows, cache.expandedCols) > 0;
        std::swap(cache.expandedRows, cache.nextExpandedRows);
        std::swap(cache.expandedCols, cache.nextExpandedCols);
        if (changed)
        {
          sortPattern(m_execution, m_cols, {}, {0}, cache.expandedRows,
            cache.expandedCols, cache.expandedOrder);
          cache.expandedStarts = cache.expandedOrder.visit([&](const auto& order) {
            return distinctStarts(
              m_execution, cache.expandedRows, cache.expandedCols, order);
          });
          cache.matrix = -1;
        }
        fillDiagonal(m_execution, entries + extra, diagonal, fixed, cache.fixed, m_owned,
          cache.expandedValues);
        if (fixed)
        {
          const size_t moved = countEliminated(m_execution, entries + extra,
            cache.expandedRows, cache.expandedCols, cache.fixed);
          eliminate(m_execution, entries + extra, cache.expandedRows, cache.expandedCols,
            cache.fixed, cache.fixedValues, cache.expandedValues, reserve(cache, moved));
        }
      }

      /**
       * @brief Reserves contributions to the right-hand side at the end of
       * the list.
       * @param count Number of contributions.
       * @returns Buffers of the list and position of the first reserved
       * contribution.
       */
      Contributions reserve(CacheType& cache, size_t count) const
      {
        Contributions res;
        res.offset = cache.contributions;
        const size_t size = res.offset + count;
        const size_t capacity = cache.contributionRows.extent(0);
        if (size > capacity)
        {
          // The list keeps its listed contributions and grows geometrically,
          // the cache keeping it for the next assemblies.
          const size_t grown = std::max(size, 2 * capacity);
          const auto alloc =
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing);
          ::Kokkos::resize(alloc, cache.contributionRows, grown);
          ::Kokkos::resize(alloc, cache.contributionValues, grown);
        }
        cache.contributions = size;
        res.rows = cache.contributionRows;
        res.values = cache.contributionValues;
        return res;
      }

      /**
       * @brief Adds the listed contributions to the right-hand side.
       *
       * The contributions to a row are summed in list order, so the
       * right-hand side does not depend on the scheduling of the kernels.
       */
      void sumContributions(CacheType& cache) const
      {
        const size_t count = cache.contributions;
        if (count == 0)
          return;
        const auto rows =
          ::Kokkos::subview(cache.contributionRows, std::make_pair(size_t(0), count));
        // The grouping depends on the listed rows only, which do not change
        // between the assemblies of a same problem.
        if (cache.groupedRows.extent(0) != count ||
          countDifferences(m_execution, rows, cache.groupedRows) > 0)
        {
          cache.groupedRows = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin grouped contribution rows")),
            count);
          ::Kokkos::deep_copy(m_execution, cache.groupedRows, rows);
          const ::Kokkos::View<Index*, MemorySpace> keys(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin contribution keys")),
            count);
          ::Kokkos::deep_copy(m_execution, keys, rows);
          cache.groupedPositions = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin contribution positions")),
            count);
          sequence(m_execution, cache.groupedPositions);
          // Without segments, the rank of a contribution is its position.
          const typename ReferenceType::Ranks ranks;
          const auto groups =
            ReferenceType::group(m_execution, keys, cache.groupedPositions, ranks);
          cache.groupStarts = groups.starts;
          cache.groupKeys = groups.keys;
          cache.vectorSegments.clear();
        }
        const typename ReferenceType::Groups groups{cache.groupStarts, cache.groupKeys};
        ReferenceType::accumulate(m_execution, cache.groupedPositions, groups,
          cache.contributionValues, cache.rhs);
      }

      /**
       * @brief Adds the entries of the linear integrators to the right-hand
       * side, as contributions listed in entry order.
       * @param offsets Position of the first entry of every integrator,
       * followed by the entry count.
       *
       * The grouping of the entries by row is kept while the integrators
       * keep their layouts and entry counts, which fix the entry indices.
       * Integrators restricted to attributes group their entries at every
       * assembly, since the attributes of the cells may change.
       */
      void sumEntries(CacheType& cache, const std::vector<size_t>& offsets) const
      {
        assert(cache.contributions == 0);
        const size_t count = offsets.back();
        if (count == 0)
          return;
        std::vector<COOSegment> segments;
        bool restricted = false;
        for (size_t k = 0; k < m_vectors.size(); ++k)
        {
          const auto& v = m_vectors[k];
          COOSegment key{v.integrator->getLayout(), v.row, 0, offsets[k + 1] - offsets[k],
            v.integrator->getSpacing(), false};
          restricted = restricted || !key.layout.attributes.empty();
          segments.push_back(std::move(key));
        }
        if (restricted || cache.vectorSegments != segments)
        {
          const ::Kokkos::View<Index*, MemorySpace> keys(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin contribution keys")),
            count);
          for (size_t k = 0; k < m_vectors.size(); ++k)
          {
            shiftEntries(m_execution, cache.entryIndices, offsets[k], offsets[k + 1],
              m_vectors[k].row, keys);
          }
          cache.groupedPositions = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin contribution positions")),
            count);
          sequence(m_execution, cache.groupedPositions);
          const typename ReferenceType::Ranks ranks;
          const auto groups =
            ReferenceType::group(m_execution, keys, cache.groupedPositions, ranks);
          cache.groupStarts = groups.starts;
          cache.groupKeys = groups.keys;
          // The listed contributions are grouped again.
          cache.groupedRows = ::Kokkos::View<Index*, MemorySpace>();
          cache.vectorSegments =
            restricted ? std::vector<COOSegment>{} : std::move(segments);
        }
        const typename ReferenceType::Groups groups{cache.groupStarts, cache.groupKeys};
        ReferenceType::accumulate(
          m_execution, cache.groupedPositions, groups, cache.entryValues, cache.rhs);
      }

      /**
       * @brief Lists buffered entries of a vector as contributions to the
       * right-hand side.
       * @param begin Position of the first entry in the buffers.
       * @param count Entry count. @param row Offset of the space.
       * @param fixed Whether degrees of freedom are fixed.
       * @param x Expansions on the masters. @param expanded Whether they apply.
       */
      void listEntries(CacheType& cache, size_t begin, size_t count, Index row,
        bool fixed, const Expansions& x, bool expanded) const
      {
        const size_t n = countVector(m_execution, begin, count, cache.entryIndices,
          cache.entryValues, row, fixed, cache.fixed, x, expanded);
        if (n == 0)
          return;
        listVector(m_execution, begin, count, cache.entryIndices, cache.entryValues, row,
          fixed, cache.fixed, x, expanded, reserve(cache, n));
      }

      /**
       * @brief Lists the entries of the linear integrators and of the
       * preassembled vectors as contributions to the right-hand side.
       */
      void integrateVectors(CacheType& cache, bool fixed, const Expansions* x)
      {
        std::vector<size_t> offsets(m_vectors.size() + 1, 0);
        for (size_t k = 0; k < m_vectors.size(); ++k)
          offsets[k + 1] = offsets[k] + m_vectors[k].integrator->prepare();
        resize(cache.entryIndices, offsets.back(), "Rodin entry indices");
        resize(cache.entryValues, offsets.back(), "Rodin entry values");
        const Expansions none;
        const Expansions& expansions = x ? *x : none;
        // Without constraints nor preassembled vectors, the entries are the
        // contributions, read in place.
        const bool direct = !fixed && !x && m_preassembled.empty();
        ::Rodin::Kokkos::Internal::clearLaunches();
        for (size_t k = 0; k < m_vectors.size(); ++k)
        {
          const auto& v = m_vectors[k];
          v.integrator->fill(
            cache.entryValues, offsets[k], v.sign, false, cache.entryIndices);
          ::Rodin::Kokkos::Internal::checkLaunches("A device linear integrator");
          if (!direct)
          {
            listEntries(cache, offsets[k], offsets[k + 1] - offsets[k], v.row, fixed,
              expansions, x != nullptr);
          }
        }
        if (direct)
        {
          sumEntries(cache, offsets);
          return;
        }
        // The integrators may still read the buffers the vectors reuse.
        if (!m_preassembled.empty())
          m_execution.fence("Rodin linear integrators");
        std::unordered_map<PetscInt, Index> locals;
        if (m_distribution && !m_preassembled.empty())
          locals = m_distribution->getLocals();
        for (const auto& p : m_preassembled)
        {
          PetscInt n = 0;
          PetscErrorCode ierr = VecGetLocalSize(p.vec, &n);
          assert(ierr == PETSC_SUCCESS);
          const size_t count = static_cast<size_t>(n);
          resize(cache.entryIndices, count, "Rodin entry indices");
          resize(cache.entryValues, count, "Rodin entry values");
          getVector(p.vec, cache.entryValues);
          Index row = p.row;
          // Alive until the fence that ends the iteration.
          std::vector<Index> indices;
          if (m_distribution)
          {
            // The owned entries are numbered locally, each once.
            PetscInt begin = 0, end = 0;
            ierr = VecGetOwnershipRange(p.vec, &begin, &end);
            assert(ierr == PETSC_SUCCESS);
            const PetscInt offset = m_distribution->getOffset(p.row);
            indices.resize(count);
            for (size_t i = 0; i < count; ++i)
            {
              const auto it = locals.find(offset + begin + static_cast<PetscInt>(i));
              assert(it != locals.end());
              indices[i] = it->second;
            }
            copy(cache.entryIndices, indices, "Rodin entry indices");
            row = 0;
          }
          else
          {
            sequence(m_execution, cache.entryIndices);
          }
          listEntries(cache, 0, count, row, fixed, expansions, x != nullptr);
          // The buffers are reused by the next vector.
          m_execution.fence("Rodin preassembled vector");
          (void)ierr;
        }
      }

      /**
       * @brief Copies a vector, wherever it stores its values.
       */
      void getVector(::Vec vec, const ::Kokkos::View<Real*, MemorySpace>& res) const
      {
        using Unmanaged = ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>;
        const size_t count = res.extent(0);
        const PetscScalar* array = nullptr;
        PetscMemType type;
        PetscErrorCode ierr = VecGetArrayReadAndMemType(vec, &array, &type);
        assert(ierr == PETSC_SUCCESS);
        if (PetscMemTypeHost(type))
        {
          ::Kokkos::deep_copy(m_execution, res,
            ::Kokkos::View<const PetscScalar*, ::Kokkos::HostSpace, Unmanaged>(
              array, count));
        }
        else
        {
          ::Kokkos::deep_copy(m_execution, res,
            ::Kokkos::View<const PetscScalar*,
              ::Kokkos::DefaultExecutionSpace::memory_space, Unmanaged>(array, count));
        }
        // PETSc may release the array once it is restored.
        m_execution.fence("Rodin preassembled vector");
        ierr = VecRestoreArrayReadAndMemType(vec, &array);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      /**
       * @brief Hands the coordinate list to the matrix.
       *
       * The matrix is laid out as a PETSc Kokkos matrix unless it already has
       * a type. The pattern is analysed again only when it is not the one the
       * matrix was last preallocated with.
       * @param expanded Whether the coordinates are expanded on the masters.
       */
      template <class Indices, class Scalars>
      void setMatrix(CacheType& cache, ::Mat A, bool expanded, const Indices& coordRows,
        const Indices& coordCols, const COOOrder<MemorySpace>& order,
        const Indices& starts, const Scalars& values) const
      {
        PetscErrorCode ierr;
        PetscInt64 id = -1;
        ierr = PetscObjectGetId(reinterpret_cast<PetscObject>(A), &id);
        assert(ierr == PETSC_SUCCESS);

        if (cache.matrix != id || cache.expanded != expanded)
        {
          PetscBool assembled = PETSC_FALSE;
          ierr = MatAssembled(A, &assembled);
          assert(ierr == PETSC_SUCCESS);
          if (!assembled)
          {
            if (m_distribution)
            {
              const auto& columns = getColumns();
              ierr = MatSetSizes(A, m_distribution->getOwnedSize(),
                columns.getOwnedSize(), m_distribution->size, columns.size);
            }
            else
            {
              const PetscInt m = static_cast<PetscInt>(m_rows);
              const PetscInt n = static_cast<PetscInt>(m_cols);
              ierr = MatSetSizes(A, m, n, m, n);
            }
            assert(ierr == PETSC_SUCCESS);
            MatType type = nullptr;
            ierr = MatGetType(A, &type);
            assert(ierr == PETSC_SUCCESS);
            if (!type)
            {
              ierr = MatSetType(A, MATAIJKOKKOS);
              assert(ierr == PETSC_SUCCESS);
            }
          }
          const size_t distinct = starts.extent(0) - 1;
          const ::Kokkos::View<PetscInt*, MemorySpace> rows(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin PETSc rows")),
            distinct);
          const ::Kokkos::View<PetscInt*, MemorySpace> cols(
            ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin PETSc cols")),
            distinct);
          order.visit([&](const auto& o) {
            listDistinct(m_execution, coordRows, coordCols, o, starts, rows, cols);
          });
          // PETSc may modify the coordinates it is given, which are temporaries.
          const auto hostRows =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), rows);
          const auto hostCols =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), cols);
          if (m_distribution)
          {
            const auto& rowGlobals = m_distribution->globals;
            const auto& colGlobals = getColumns().globals;
            for (size_t i = 0; i < distinct; ++i)
            {
              hostRows(i) = rowGlobals[hostRows(i)];
              hostCols(i) = colGlobals[hostCols(i)];
            }
          }
          ierr = MatSetPreallocationCOO(
            A, static_cast<PetscCount>(distinct), hostRows.data(), hostCols.data());
          assert(ierr == PETSC_SUCCESS);
          cache.matrix = id;
          cache.expanded = expanded;
        }

        resize(cache.merged, starts.extent(0) - 1, "Rodin COO merged values");
        order.visit(
          [&](const auto& o) { merge(m_execution, o, starts, values, cache.merged); });
        m_execution.fence("Rodin COO merge");

        PetscBool device = PETSC_FALSE;
        ierr = PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(A), &device,
          MATSEQAIJKOKKOS, MATMPIAIJKOKKOS, MATSEQAIJCUSPARSE, MATMPIAIJCUSPARSE, "");
        assert(ierr == PETSC_SUCCESS);
        if (device ||
          ::Kokkos::SpaceAccessibility<::Kokkos::HostSpace, MemorySpace>::accessible)
        {
          ierr = MatSetValuesCOO(A, cache.merged.data(), INSERT_VALUES);
        }
        else
        {
          const auto host =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), cache.merged);
          ierr = MatSetValuesCOO(A, host.data(), INSERT_VALUES);
        }
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      /**
       * @brief Copies the right-hand side where the vector stores it.
       */
      void setVector(::Vec b, const ::Kokkos::View<Real*, MemorySpace>& rhs) const
      {
        using Unmanaged = ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>;
        const auto source = ::Kokkos::subview(rhs, std::make_pair(size_t(0), m_rows));
        PetscScalar* array = nullptr;
        PetscMemType type;
        PetscErrorCode ierr = VecGetArrayWriteAndMemType(b, &array, &type);
        assert(ierr == PETSC_SUCCESS);
        if (PetscMemTypeHost(type))
        {
          ::Kokkos::deep_copy(
            ::Kokkos::View<PetscScalar*, ::Kokkos::HostSpace, Unmanaged>(array, m_rows),
            source);
        }
        else
        {
          // Ordered after the kernels that sum the right-hand side, waited
          // for before finish() returns.
          ::Kokkos::deep_copy(m_execution,
            ::Kokkos::View<PetscScalar*, ::Kokkos::DefaultExecutionSpace::memory_space,
              Unmanaged>(array, m_rows),
            source);
        }
        ierr = VecRestoreArrayWriteAndMemType(b, &array);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      /// @brief Gets the distribution of the columns. @return Distribution of the columns.
      const COODistribution& getColumns() const
      {
        assert(m_distribution);
        return m_columns ? *m_columns : *m_distribution;
      }

      /**
       * @brief Whether a process fixes a degree of freedom; collective.
       * @param local Whether this process fixes one.
       */
      bool isFixed(bool local) const
      {
        int in = local ? 1 : 0, out = 0;
        MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_LOR, m_distribution->comm);
        return out != 0;
      }

      /**
       * @brief Sums a value per local degree of freedom over the processes
       * holding it; collective.
       * @param[in,out] local Value per local degree of freedom, replaced by
       * the sum over the processes.
       */
      void sum(std::vector<PetscScalar>& local) const
      {
        const auto& dist = *m_distribution;
        const PetscInt owned = dist.getOwnedSize();
        std::vector<PetscInt> ghosts;
        std::vector<PetscInt> position(m_rows);
        for (size_t i = 0; i < m_rows; ++i)
        {
          const PetscInt g = dist.globals[i];
          if (g >= dist.begin && g < dist.end)
          {
            position[i] = g - dist.begin;
          }
          else
          {
            position[i] = owned + static_cast<PetscInt>(ghosts.size());
            ghosts.push_back(g);
          }
        }
        ::Vec v;
        PetscErrorCode ierr = VecCreateGhost(dist.comm, owned, dist.size,
          static_cast<PetscInt>(ghosts.size()), ghosts.data(), &v);
        assert(ierr == PETSC_SUCCESS);
        ::Vec form;
        ierr = VecGhostGetLocalForm(v, &form);
        assert(ierr == PETSC_SUCCESS);
        PetscScalar* array = nullptr;
        ierr = VecGetArray(form, &array);
        assert(ierr == PETSC_SUCCESS);
        for (size_t i = 0; i < m_rows; ++i)
          array[position[i]] = local[i];
        ierr = VecRestoreArray(form, &array);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateBegin(v, ADD_VALUES, SCATTER_REVERSE);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateEnd(v, ADD_VALUES, SCATTER_REVERSE);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGetArray(form, &array);
        assert(ierr == PETSC_SUCCESS);
        for (size_t i = 0; i < m_rows; ++i)
          local[i] = array[position[i]];
        ierr = VecRestoreArray(form, &array);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostRestoreLocalForm(v, &form);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDestroy(&v);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      /**
       * @brief Makes the fixed degrees of freedom the same on every process
       * holding them; collective.
       *
       * A degree of freedom is fixed when one process fixes it, to the mean
       * of the values the processes fixing it impose.
       * @param[in,out] mask Whether each local degree of freedom is fixed.
       * @param[in,out] values Fixed value of each local degree of freedom.
       */
      void synchronize(std::vector<unsigned char>& mask, std::vector<Real>& values) const
      {
        std::vector<PetscScalar> count(m_rows), value(m_rows);
        for (size_t i = 0; i < m_rows; ++i)
        {
          count[i] = mask[i] ? PetscScalar(1) : PetscScalar(0);
          value[i] = mask[i] ? PetscScalar(values[i]) : PetscScalar(0);
        }
        sum(count);
        sum(value);
        for (size_t i = 0; i < m_rows; ++i)
        {
          const Real c = PetscRealPart(count[i]);
          mask[i] = c > Real(0.5);
          values[i] = mask[i] ? PetscRealPart(value[i]) / c : Real(0);
        }
      }

      /**
       * @brief Resolves the identifications of a distributed system;
       * collective.
       *
       * The fixed degrees of freedom are first synchronized, then every
       * condition is resolved in a map of the global degrees of freedom.
       * The masters the process does not hold become local degrees of
       * freedom past the held ones, with their global index, so that the
       * entries the expansion moves to them are sent to their owners.
       */
      typename ConstraintsType::Resolved resolveDistributed()
      {
#ifndef RODIN_USE_MPI
        Alert::Exception() << "Distributed device systems require MPI." << Alert::Raise;
        return {};
#else
        using ConstraintMap = ::Rodin::PETSc::Assembly::ConstraintMap;
        auto& dist = *m_distribution;
        assert(m_rows == m_cols && !m_columns);
        auto mask = m_constraints.getFixed();
        auto values = m_constraints.getValues();
        mask.resize(m_rows, 0);
        values.resize(m_rows, Real(0));
        if (isFixed(m_constraints.isConstrained()))
          synchronize(mask, values);
        m_constraints.setFixed(std::move(mask), std::move(values));
        std::vector<Index> globals(dist.globals.begin(), dist.globals.end());
        ConstraintMap map(static_cast<size_t>(dist.size));
        auto res = m_constraints.resolve(
          m_execution, map, globals, dist.owned,
          [&](size_t start) {
            return static_cast<Index>(dist.getOffset(static_cast<Index>(start)));
          },
          [&](ConstraintMap& m) { m.finalize(dist.comm); });
        dist.globals.assign(globals.begin(), globals.end());
        m_rows = m_cols = globals.size();
        copy(m_owned, dist.owned, "Rodin owned DOFs");
        return res;
#endif
      }

      /**
       * @brief Sends the right-hand side to the owners of its rows.
       *
       * The vector is preallocated with the global index of each local row,
       * once per vector, and PETSc sums the entries of the rows held by
       * several processes.
       */
      void setDistributedVector(CacheType& cache, ::Vec b) const
      {
        PetscErrorCode ierr;
        PetscInt64 id = -1;
        ierr = PetscObjectGetId(reinterpret_cast<PetscObject>(b), &id);
        assert(ierr == PETSC_SUCCESS);
        if (cache.vector != id)
        {
          // PETSc may modify the coordinates it is given.
          std::vector<PetscInt> globals = m_distribution->globals;
          ierr =
            VecSetPreallocationCOO(b, static_cast<PetscCount>(m_rows), globals.data());
          assert(ierr == PETSC_SUCCESS);
          cache.vector = id;
        }
        const auto source =
          ::Kokkos::subview(cache.rhs, std::make_pair(size_t(0), m_rows));
        PetscBool device = PETSC_FALSE;
        ierr = PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(b), &device,
          VECSEQKOKKOS, VECMPIKOKKOS, VECSEQCUDA, VECMPICUDA, "");
        assert(ierr == PETSC_SUCCESS);
        if (device ||
          ::Kokkos::SpaceAccessibility<::Kokkos::HostSpace, MemorySpace>::accessible)
        {
          ierr = VecSetValuesCOO(b, source.data(), INSERT_VALUES);
        }
        else
        {
          const auto host =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), source);
          ierr = VecSetValuesCOO(b, host.data(), INSERT_VALUES);
        }
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      ExecutionSpace m_execution;
      size_t m_rows;
      size_t m_cols;
      ConstraintsType m_constraints;
      std::vector<Matrix> m_matrices;
      std::vector<Vector> m_vectors;
      std::vector<Operator> m_operators;
      std::vector<Preassembled> m_preassembled;
      std::optional<COODistribution> m_distribution;
      std::optional<COODistribution> m_columns;
      ::Kokkos::View<unsigned char*, MemorySpace> m_owned;
  };

  /**
   * @brief Prepares the PETSc objects of a distributed device assembly.
   * @param axb PETSc linear system.
   * @param rows Number of owned rows. @param globalRows Global number of rows.
   * @param cols Number of owned columns. @param globalCols Global number of columns.
   * @param doVector Whether the right-hand side is assembled.
   */
  inline void prepareDeviceSystem(PETSc::Math::LinearSystem& axb, size_t rows,
    size_t globalRows, size_t cols, size_t globalCols, bool doVector)
  {
    PetscErrorCode ierr;
    const PetscInt m = static_cast<PetscInt>(rows);
    const PetscInt n = static_cast<PetscInt>(cols);
    const PetscInt M = static_cast<PetscInt>(globalRows);
    const PetscInt N = static_cast<PetscInt>(globalCols);
    if (doVector)
    {
      ierr = PETSc::Assembly::VectorSetup(axb.getVector())
               .prepare({m, M, VECKOKKOS, false, false});
      assert(ierr == PETSC_SUCCESS);
    }
    // The solution carries the previous iterate as the initial guess.
    ierr = PETSc::Assembly::VectorSetup(axb.getSolution())
             .prepare({n, N, VECKOKKOS, false, false});
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  /**
   * @brief Prepares the PETSc objects of a device assembly.
   * @param axb PETSc linear system.
   * @param rows Number of rows. @param cols Number of columns.
   * @param doVector Whether the right-hand side is assembled.
   */
  inline void prepareDeviceSystem(
    PETSc::Math::LinearSystem& axb, size_t rows, size_t cols, bool doVector)
  {
    prepareDeviceSystem(axb, rows, rows, cols, cols, doVector);
  }
}

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a PETSc bilinear form on a Kokkos mesh.
   *
   * The integrators of the form, possibly negated, are integrated on the
   * device and handed to PETSc as a coordinate list in device memory, as in
   * the problem assemblies. The pattern holds the diagonal.
   *
   * @tparam Solution Solution variable type.
   * @tparam TrialFES Trial space on a resident Kokkos mesh.
   * @tparam TestFES Test space on the same resident Kokkos mesh.
   */
  template <class Solution, class TrialFES, class TestFES>
  class Kokkos<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat,
        Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
  {
    public:
      /// @brief Operator type.
      using OperatorType = ::Mat;

      /// @brief Bilinear form type.
      using BilinearFormType =
        Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the trial mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<TrialFES>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::DeviceCOOSystem<BackendType>;

      /// @brief Default constructor.
      Kokkos() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Kokkos(const Kokkos& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Kokkos(Kokkos&& other)
        : Parent(std::move(other)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the form on the device and fills its matrix.
       * @param res Matrix of the form.
       * @param input Input data.
       */
      void execute(OperatorType& res, const InputType& input) override
      {
        const auto& trial = input.getTrialFES();
        const auto& test = input.getTestFES();
        constexpr bool distributed = Context::IsMPI<typename FormLanguage::Traits<
          typename FormLanguage::Traits<TrialFES>::MeshType>::ContextType>::value;
        // The kernels number the local degrees of freedom, ghosts included.
        size_t rows, cols;
        if constexpr (distributed)
        {
          rows = test.getLocalSize();
          cols = trial.getLocalSize();
        }
        else
        {
          rows = test.getSize();
          cols = trial.getSize();
        }
        SystemType system(trial.getMesh().getContext().getExecutionSpace(), rows, cols);
        if constexpr (distributed)
        {
          system.distribute(Internal::COODistribution::of(test),
            Internal::COODistribution::of(trial));
        }
        for (const auto& lbfi : input.getLocalBFIs())
          system.addMatrix(lbfi, 0, 0);
        for (const auto& gbfi : input.getGlobalBFIs())
          system.addMatrix(gbfi, 0, 0);
        system.finish(m_cache, res, nullptr);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }

    private:
      typename SystemType::CacheType m_cache;
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a PETSc linear form on a Kokkos mesh.
   *
   * The integrators of the form, possibly negated, are integrated on the
   * device and summed into a PETSc Kokkos vector.
   *
   * @tparam FES Test space on a resident Kokkos mesh.
   */
  template <class FES>
  class Kokkos<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
  {
    public:
      /// @brief Vector type.
      using VectorType = ::Vec;

      /// @brief Linear form type.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<VectorType, LinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the test mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<FES>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::DeviceCOOSystem<BackendType>;

      /// @brief Default constructor.
      Kokkos() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Kokkos(const Kokkos& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Kokkos(Kokkos&& other)
        : Parent(std::move(other)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached buffers.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the form on the device and fills its vector.
       * @param res Vector of the form.
       * @param input Input data.
       */
      void execute(VectorType& res, const InputType& input) override
      {
        const auto& fes = input.getFES();
        constexpr bool distributed = Context::IsMPI<typename FormLanguage::Traits<
          typename FormLanguage::Traits<FES>::MeshType>::ContextType>::value;
        PetscInt n = static_cast<PetscInt>(fes.getSize());
        const PetscInt N = n;
        size_t size = fes.getSize();
        if constexpr (distributed)
        {
          // The kernels number the local degrees of freedom, ghosts included.
          size = fes.getLocalSize();
          Index begin = 0, end = 0;
          fes.getOwnershipRange(begin, end);
          n = static_cast<PetscInt>(end - begin);
        }
        SystemType system(fes.getMesh().getContext().getExecutionSpace(), size, 0);
        if constexpr (distributed)
          system.distribute(Internal::COODistribution::of(fes));
        for (const auto& lfi : input.getLFIs())
          system.addVector(lfi, 0, Real(1));
        PetscErrorCode ierr =
          PETSc::Assembly::VectorSetup(res).prepare({n, N, VECKOKKOS, false, false});
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        system.finish(m_cache, nullptr, res);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }

    private:
      typename SystemType::CacheType m_cache;
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a single field PETSc problem on a Kokkos mesh.
   *
   * The integrators of the problem, possibly negated, are integrated on the
   * device and handed to PETSc as coordinate lists in device memory. The
   * Dirichlet and periodic conditions are imposed on the device, as the
   * host assemblies do.
   *
   * The entries of the preassembled forms are copied from their PETSc
   * objects and combined with the integrated ones.
   *
   * @tparam U Trial function on a resident Kokkos mesh.
   * @tparam V Test function on the same mesh.
   */
  template <class U, class V>
  class Kokkos<PETSc::Math::LinearSystem,
    Variational::Problem<PETSc::Math::LinearSystem, U, V>>
    final : public AssemblyBase<PETSc::Math::LinearSystem,
              Variational::Problem<PETSc::Math::LinearSystem, U, V>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U, V>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the trial mesh context.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U>::FESType>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::DeviceCOOSystem<BackendType>;

      /// @brief Default constructor.
      Kokkos() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Kokkos(const Kokkos& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Kokkos(Kokkos&& other)
        : Parent(std::move(other)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the problem on the device and builds its system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       */
      void execute(LinearSystemType& axb, const InputType& input) override
      {
        assemble(axb, input, true, true);
      }

      /**
       * @brief Assembles one side of the system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        const bool lhs = target == Variational::AssemblyTarget::LHS;
        assemble(axb, input, lhs, !lhs);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }

    private:
      /// @brief Whether the mesh is distributed over MPI processes.
      static constexpr bool IsDistributed =
        Context::IsMPI<typename FormLanguage::Traits<typename FormLanguage::Traits<
          typename FormLanguage::Traits<U>::FESType>::MeshType>::ContextType>::value;

      void assemble(
        LinearSystemType& axb, const InputType& input, bool doMatrix, bool doVector)
      {
        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        const auto& trial = u.getFiniteElementSpace();
        const auto& test = v.getFiniteElementSpace();
        // The kernels number the local degrees of freedom, ghosts included.
        size_t rows, cols;
        if constexpr (IsDistributed)
        {
          rows = test.getLocalSize();
          cols = trial.getLocalSize();
        }
        else
        {
          rows = test.getSize();
          cols = trial.getSize();
        }
        SystemType system(trial.getMesh().getContext().getExecutionSpace(), rows, cols);
        if constexpr (IsDistributed)
        {
          if (&trial.getMesh() != &test.getMesh() || rows != cols)
          {
            Alert::Exception()
              << "A distributed device problem requires the trial and test "
                 "functions to share their space."
              << Alert::Raise;
          }
          system.distribute(Internal::COODistribution::of(trial));
        }
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() == u.getUUID())
            system.fix(dbc, 0);
        }
        for (auto& pbc : pb.getPBCs())
        {
          if (pbc.getOperand().getUUID() == u.getUUID())
            system.identify(pbc, 0);
        }
        // The matrix is still needed on the right-hand side to lift the
        // constrained values.
        if (doMatrix || system.hasConstraints())
        {
          for (const auto& lbfi : pb.getLocalBFIs())
            system.addMatrix(lbfi, 0, 0);
          for (const auto& gbfi : pb.getGlobalBFIs())
            system.addMatrix(gbfi, 0, 0);
          for (const auto& bf : pb.getBFs())
            system.addMatrix(bf.getOperator(), 0, 0);
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        if (doVector)
        {
          for (const auto& lfi : pb.getLFIs())
            system.addVector(lfi, 0, Real(-1));
          for (const auto& lf : pb.getLFs())
            system.addVector(lf.getVector(), 0);
        }
        if constexpr (IsDistributed)
        {
          Index begin = 0, end = 0;
          trial.getOwnershipRange(begin, end);
          const size_t owned = end - begin;
          Internal::prepareDeviceSystem(
            axb, owned, test.getSize(), owned, trial.getSize(), doVector);
        }
        else
        {
          Internal::prepareDeviceSystem(axb, rows, cols, doVector);
        }
        system.finish(m_cache, doMatrix ? axb.getOperator() : nullptr,
          doVector ? axb.getVector() : nullptr);
      }

      typename SystemType::CacheType m_cache;
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a mixed PETSc problem on a Kokkos mesh.
   *
   * Each field occupies a block of the system, at the offsets given by the
   * problem input, and the system is assembled as in the single field case.
   */
  template <class U1, class U2, class U3, class... Us>
  class Kokkos<PETSc::Math::LinearSystem,
    Variational::Problem<PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
    final : public AssemblyBase<PETSc::Math::LinearSystem,
              Variational::Problem<PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the mesh context of the first field.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::DeviceCOOSystem<BackendType>;

      /// @brief Default constructor.
      Kokkos() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Kokkos(const Kokkos& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Kokkos(Kokkos&& other)
        : Parent(std::move(other)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the problem on the device and builds its system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       */
      void execute(LinearSystemType& axb, const InputType& input) override
      {
        assemble(axb, input, true, true);
      }

      /**
       * @brief Assembles one side of the system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        const bool lhs = target == Variational::AssemblyTarget::LHS;
        assemble(axb, input, lhs, !lhs);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }

    private:
      /// @brief Whether the mesh of the first field is distributed over MPI processes.
      static constexpr bool IsDistributed =
        Context::IsMPI<typename FormLanguage::Traits<typename FormLanguage::Traits<
          typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType>::value;

      void assemble(
        LinearSystemType& axb, const InputType& input, bool doMatrix, bool doVector)
      {
        auto& pb = input.getProblemBody();
        const auto& us = input.getTrialFunctions();
        const auto& vs = input.getTestFunctions();
        const auto& trialUUIDMap = input.getTrialUUIDMap();
        const auto& testUUIDMap = input.getTestUUIDMap();
        // The kernels number the local degrees of freedom, ghosts included,
        // and the blocks of a distributed system follow each other in the
        // local numbering.
        std::vector<Index> trialOffsets(input.getTrialOffsets().begin(),
          input.getTrialOffsets().end());
        std::vector<Index> testOffsets(
          input.getTestOffsets().begin(), input.getTestOffsets().end());
        size_t rows = input.getTotalTestSize();
        size_t cols = input.getTotalTrialSize();
        Internal::COODistribution distribution;
        if constexpr (IsDistributed)
        {
          std::vector<size_t> sizes(trialOffsets.size());
          std::vector<size_t> globalSizes(trialOffsets.size());
          cols = 0;
          size_t i = 0;
          us.apply([&](const auto& u) {
            const auto& fes = u.get().getFiniteElementSpace();
            trialOffsets[i] = static_cast<Index>(cols);
            sizes[i] = fes.getLocalSize();
            globalSizes[i] = fes.getSize();
            cols += sizes[i];
            distribution.append(fes, static_cast<PetscInt>(input.getTrialOffsets()[i]));
            ++i;
          });
          rows = 0;
          bool same = testOffsets.size() == trialOffsets.size();
          i = 0;
          vs.apply([&](const auto& v) {
            const auto& fes = v.get().getFiniteElementSpace();
            testOffsets[i] = static_cast<Index>(rows);
            rows += fes.getLocalSize();
            same = same && i < sizes.size() && fes.getLocalSize() == sizes[i] &&
              fes.getSize() == globalSizes[i];
            ++i;
          });
          if (!same)
          {
            Alert::Exception()
              << "A distributed device problem requires each test function to "
                 "share the space of the trial function of its block."
              << Alert::Raise;
          }
          distribution.close(
            us.template get<0>().get().getFiniteElementSpace().getMesh().getContext()
              .getCommunicator(),
            static_cast<PetscInt>(input.getTotalTrialSize()));
        }
        const auto trialOffset = [&](const auto& uuid) -> Index {
          const auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return trialOffsets[it->second];
        };
        const auto testOffset = [&](const auto& uuid) -> Index {
          const auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return testOffsets[it->second];
        };

        const auto& u1 = us.template get<0>().get();
        SystemType system(
          u1.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(), rows,
          cols);
        if constexpr (IsDistributed)
          system.distribute(std::move(distribution));
        for (auto& dbc : pb.getDBCs())
        {
          const Index offset = trialOffset(dbc.getOperand().getUUID());
          const auto master = dbc.getValueUUID();
          system.fix(dbc, offset, master ? trialOffset(*master) : offset);
        }
        for (auto& pbc : pb.getPBCs())
          system.identify(pbc, trialOffset(pbc.getOperand().getUUID()));
        if (doMatrix || system.hasConstraints())
        {
          for (const auto& lbfi : pb.getLocalBFIs())
          {
            system.addMatrix(lbfi, testOffset(lbfi.getTestFunction().getUUID()),
              trialOffset(lbfi.getTrialFunction().getUUID()));
          }
          for (const auto& gbfi : pb.getGlobalBFIs())
          {
            system.addMatrix(gbfi, testOffset(gbfi.getTestFunction().getUUID()),
              trialOffset(gbfi.getTrialFunction().getUUID()));
          }
          for (const auto& bf : pb.getBFs())
          {
            system.addMatrix(bf.getOperator(), testOffset(bf.getTestFunction().getUUID()),
              trialOffset(bf.getTrialFunction().getUUID()));
          }
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        if (doVector)
        {
          for (const auto& lfi : pb.getLFIs())
            system.addVector(lfi, testOffset(lfi.getTestFunction().getUUID()), Real(-1));
          for (const auto& lf : pb.getLFs())
            system.addVector(lf.getVector(), testOffset(lf.getTestFunction().getUUID()));
        }
        if constexpr (IsDistributed)
        {
          const size_t owned = system.getDistribution().getOwnedSize();
          Internal::prepareDeviceSystem(axb, owned, input.getTotalTestSize(), owned,
            input.getTotalTrialSize(), doVector);
        }
        else
        {
          Internal::prepareDeviceSystem(axb, rows, cols, doVector);
        }
        system.finish(m_cache, doMatrix ? axb.getOperator() : nullptr,
          doVector ? axb.getVector() : nullptr);
      }

      typename SystemType::CacheType m_cache;
  };
}

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Gathers the local form of a ghosted PETSc vector in the local
   * numbering of the device degrees of freedom.
   *
   * Kernel functor at namespace scope, as NVCC requires.
   */
  template <class MemorySpace>
  struct GhostedCoefficientsGather
  {
    /// @brief Coefficients in the local numbering.
    ::Kokkos::View<Real*, MemorySpace> res;

    /// @brief Local form of the vector, owned entries then ghosts.
    ::Kokkos::View<const Real*, MemorySpace> form;

    /// @brief Position in the local form of each local degree of freedom.
    ::Kokkos::View<const Index*, MemorySpace> position;

    /// @param i Local degree of freedom.
    KOKKOS_INLINE_FUNCTION void operator()(Index i) const
    {
      res(i) = form(position(i));
    }
  };

  /**
   * @brief Copies the coefficients of a PETSc grid function to the device.
   *
   * The values are read where PETSc keeps them, on the host or on the
   * device. On a distributed mesh the ghosted local form is read, without
   * communication: the collective writes of PETSc grid functions keep the
   * ghosts up to date. Its owned entries are numbered from the first owned
   * global index and its ghosts follow in the order of the local degrees of
   * freedom, as PETSc::Variational::GridFunction lays them out.
   */
  template <>
  struct GridFunctionCoefficients<::Vec>
  {
    /**
     * @param fes Device space of the grid function.
     * @param data Coefficients of the grid function.
     * @returns Device coefficients in the local numbering, after transfer
     * completion.
     */
    template <class FES>
    static ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> upload(
      const FES& fes, ::Vec data)
    {
      using MemorySpace = typename FES::MeshType::MemorySpace;
      using Unmanaged = ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>;
      using ContextType = typename FES::MeshType::ContextType;
      const auto& execution = fes.getMesh().getContext().getExecutionSpace();

      ::Vec form = data;
      PetscErrorCode ierr = PETSC_SUCCESS;
      if constexpr (::Rodin::Context::IsMPI<ContextType>::value)
      {
        ierr = VecGhostGetLocalForm(data, &form);
        assert(ierr == PETSC_SUCCESS);
        if (!form)
        {
          Alert::Exception()
            << "The PETSc vector of a distributed grid function has no ghosts."
            << Alert::Raise;
        }
      }
      PetscInt count = 0;
      ierr = VecGetLocalSize(form, &count);
      assert(ierr == PETSC_SUCCESS);

      ::Kokkos::View<Real*, MemorySpace> values(
        ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
          std::string("Rodin PETSc grid function")),
        count);
      const PetscScalar* array = nullptr;
      PetscMemType type;
      ierr = VecGetArrayReadAndMemType(form, &array, &type);
      assert(ierr == PETSC_SUCCESS);
      if (PetscMemTypeHost(type))
      {
        ::Kokkos::deep_copy(execution, values,
          ::Kokkos::View<const PetscScalar*, ::Kokkos::HostSpace, Unmanaged>(array, count));
      }
      else
      {
        ::Kokkos::deep_copy(execution, values,
          ::Kokkos::View<const PetscScalar*,
            ::Kokkos::DefaultExecutionSpace::memory_space, Unmanaged>(array, count));
      }
      // PETSc may release the array once it is restored.
      execution.fence("Rodin PETSc grid function transfer");
      ierr = VecRestoreArrayReadAndMemType(form, &array);
      assert(ierr == PETSC_SUCCESS);

      if constexpr (::Rodin::Context::IsMPI<ContextType>::value)
      {
        ierr = VecGhostRestoreLocalForm(data, &form);
        assert(ierr == PETSC_SUCCESS);

        Index begin = 0, end = 0;
        fes.getOwnershipRange(begin, end);
        const Index owned = end - begin;
        const size_t size = fes.getLocalSize();
        Index ghosts = 0;
        auto position = ::Rodin::Kokkos::Internal::upload<Index>(
          execution, size, "Rodin PETSc grid function positions", [&](Index* out) {
            for (size_t i = 0; i < size; ++i)
            {
              const Index g = fes.getGlobalIndex(static_cast<Index>(i));
              out[i] = (g >= begin && g < end) ? g - begin : owned + ghosts++;
            }
          });
        if (owned + ghosts != static_cast<Index>(count))
        {
          Alert::Exception()
            << "The PETSc vector of a distributed grid function does not have "
            << "the layout of its space." << Alert::Raise;
        }
        ::Kokkos::View<Real*, MemorySpace> res(
          ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
            std::string("Rodin grid function")),
          size);
        ::Kokkos::parallel_for("Rodin PETSc grid function gather",
          ::Kokkos::RangePolicy<std::decay_t<decltype(execution)>>(execution, 0, size),
          GhostedCoefficientsGather<MemorySpace>{ res, values, position });
        execution.fence("Rodin PETSc grid function gather");
        (void)ierr;
        return res;
      }
      else
      {
        (void)ierr;
        return values;
      }
    }
  };
}

#endif
