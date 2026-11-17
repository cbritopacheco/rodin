/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_ASSEMBLY_COOASSEMBLY_H
#define RODIN_PETSC_KOKKOS_ASSEMBLY_COOASSEMBLY_H

/**
 * @file
 * @brief Device COO assembly and distributed constraint integration for PETSc.
 */

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>

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

#include "COOCache.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Device system assembled into a coordinate list for PETSc.
   *
   * The integrators write their entries into a coordinate list in device
   * memory. The fixed degrees of freedom are eliminated and the identified
   * ones are expanded on their masters on the device, and the values summed
   * per distinct coordinate are handed to MatSetValuesCOO. The system may
   * be distributed over MPI processes.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class COOAssembly
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space of the kernels.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /**
       * @brief Range of the kernels, indexed by Index.
       *
       * The default index of a Kokkos 5 range is the index_type of the execution
       * space, a 32 bit int on the GPUs, and a bound past it aborts. One process
       * can hold more than 2^31 entries.
       */
      using Range = ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>;

      /// @brief Memory space of the system.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Cache type.
      using CacheType = COOCache<MemorySpace>;

      /// @brief Constraints of the system.
      using ConstraintsType = Constraints<Backend>;

      /// @brief Expansions of the degrees of freedom on their masters.
      using Expansions = typename ConstraintsType::Expansions;

      /// @brief Device system whose grouping sums the right-hand side.
      using ReferenceType = SystemAssembly<Backend>;

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
      using BilinearType = ::Rodin::Kokkos::Variational::BilinearFormIntegrator<Backend>;

      /// @brief Device linear integrator.
      using LinearType = ::Rodin::Kokkos::Variational::LinearFormIntegrator<Backend>;

      static_assert(std::is_same_v<PetscScalar, Real>,
        "Device assembly of PETSc problems requires real PETSc scalars.");

      /**
       * @brief Creates an empty system.
       * @param execution Execution space instance of the mesh.
       * @param rows Number of rows, the size of the test spaces.
       * @param cols Number of columns, the size of the trial spaces.
       */
      COOAssembly(const ExecutionSpace& execution, size_t rows, size_t cols)
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
       * @return Whether a degree of freedom is fixed or identified, on any
       * process once the system is distributed.
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
      void addMatrix(
        Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Index row, Index col)
      {
        const auto [integrator, sign] = getIntegrator<Backend>(lbfi);
        m_matrices.push_back({integrator, row, col, sign});
      }

      /**
       * @brief Adds a global integrator to the system.
       * @param gbfi Global integrator, possibly negated.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(
        Variational::GlobalBilinearFormIntegratorBase<Real>& gbfi, Index row, Index col)
      {
        const auto [integrator, sign] = getIntegrator<Backend>(gbfi);
        m_matrices.push_back({integrator, row, col, sign});
      }

      /**
       * @brief Adds a linear integrator to the system.
       * @param lfi Integrator, possibly negated.
       * @param row Offset of its test space.
       * @param sign Factor of the contributions.
       */
      void addVector(
        Variational::LinearFormIntegratorBase<Real>& lfi, Index row, Real sign)
      {
        const auto [integrator, s] = getIntegrator<Backend>(lfi, sign);
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
          constraints = resolveDistributed(cache);
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
          synchronize(cache, mask, values);
          copy(cache.fixed, mask, "Rodin fixed DOFs");
          copy(cache.fixedValues, values, "Rodin fixed values");
          keepFixed(cache, mask);
        }
        else if (fixed || expanded)
        {
          copy(cache.fixed, m_constraints.getFixed(), "Rodin fixed DOFs");
          copy(cache.fixedValues, m_constraints.getValues(), "Rodin fixed values");
          keepFixed(cache, m_constraints.getFixed());
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
          {
            restoreCoordinates(cache);
            expandMatrices(cache, constraints, fixed);
          }
          else if (m_operators.empty())
          {
            // The next assemblies keep the pattern, or list only its
            // restricted segments again.
            releaseCoordinates(cache);
          }
        }
        integrateVectors(cache, fixed, expanded ? &constraints.expansions : nullptr);
        sumContributions(cache);
        if (expanded)
        {
          imposeFixed(
            m_execution, m_rows, cache.imposed, cache.fixedValues, m_owned, cache.rhs);
        }
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
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Index row; ///< Offset of the test space, added to the rows.
        Index col; ///< Offset of the trial space, added to the columns.

        /**
         * @brief Adds the offsets to the coordinates of the entry.
         * @param e Position of the entry in the list.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          rows(e) += row;
          cols(e) += col;
        }
      };

      /**
       * @brief Adds the block offsets to the coordinates of a segment.
       * @param execution Execution space instance of the kernel.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param offset Position of the first entry of the segment.
       * @param count Number of entries of the segment.
       * @param row Offset of the test space.
       * @param col Offset of the trial space.
       */
      template <class Indices>
      static void shift(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, size_t offset, size_t count, Index row, Index col)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::shift",
          Range(execution, offset, offset + count), Shift<Indices>{rows, cols, row, col});
      }

      /// @brief Lists the diagonal entry i after the given offset.
      template <class Indices>
      struct DiagonalList
      {
        size_t offset; ///< Position of the first diagonal entry in the list.
        Indices rows; ///< Receives the row coordinates.
        Indices cols; ///< Receives the column coordinates.

        /**
         * @brief Lists the diagonal entry.
         * @param i Row and column of the diagonal entry.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          rows(offset + i) = i;
          cols(offset + i) = i;
        }
      };

      /**
       * @brief Lists every diagonal entry.
       * @param execution Execution space instance of the kernel.
       * @param offset Position of the first diagonal entry in the list.
       * @param count Number of diagonal entries.
       * @param rows Receives the row coordinates.
       * @param cols Receives the column coordinates.
       */
      template <class Indices>
      static void listDiagonal(const ExecutionSpace& execution, size_t offset,
        size_t count, const Indices& rows, const Indices& cols)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::diagonal",
          Range(execution, 0, count), DiagonalList<Indices>{offset, rows, cols});
      }

      /// @brief Writes the diagonal value of row i, one when it is fixed.
      template <class Fixed, class Owned, class Scalars>
      struct DiagonalFill
      {
        size_t offset; ///< Position of the first diagonal value in the list.
        bool hasFixed; ///< Whether the system has fixed degrees of freedom.
        Fixed fixed; ///< Flags of the fixed degrees of freedom.

        /**
         * @brief Flags of the degrees of freedom the process owns, empty when the system
         * is not distributed.
         */
        Owned owned;
        Scalars coo; ///< Values of the coordinate list.
        bool distributed; ///< Whether the system is distributed.

        /**
         * @brief Writes the diagonal value of the row.
         * @param i Row of the diagonal value.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          const bool one = hasFixed && fixed(i) && (!distributed || owned(i));
          coo(offset + i) = one ? PetscScalar(1) : PetscScalar(0);
        }
      };

      /**
       * @brief Writes the diagonal values, one on fixed rows.
       * @param execution Execution space instance of the kernel.
       * @param offset Position of the first diagonal value in the list.
       * @param count Number of diagonal values.
       * @param hasFixed Whether the system has fixed degrees of freedom.
       * @param fixed Flags of the fixed degrees of freedom.
       * @param owned Flags of the degrees of freedom the process owns, empty
       * when the system is not distributed.
       * @param coo Values of the coordinate list.
       *
       * On a distributed system, the owner of a fixed row alone writes one.
       */
      template <class Fixed, class Owned, class Scalars>
      static void fillDiagonal(const ExecutionSpace& execution, size_t offset,
        size_t count, bool hasFixed, const Fixed& fixed, const Owned& owned,
        const Scalars& coo)
      {
        const bool distributed = owned.extent(0) > 0;
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::diagonal",
          Range(execution, 0, count),
          DiagonalFill<Fixed, Owned, Scalars>{
            offset, hasFixed, fixed, owned, coo, distributed});
      }

      /// @brief Counts entry e when its value moves to the right-hand side.
      template <class Indices, class Fixed>
      struct EliminatedCount
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Fixed fixed; ///< Flags of the fixed degrees of freedom.

        /**
         * @brief Counts the entry when its value moves to the right-hand side.
         * @param e Position of the entry in the list.
         * @param sum Number of moved values, incremented.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += !fixed(rows(e)) && fixed(cols(e));
        }
      };

      /**
       * @brief Counts the matrix values moved to the right-hand side by the
       * columns of the fixed degrees of freedom.
       * @param execution Execution space instance of the kernel.
       * @param count Number of entries in the list.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param fixed Flags of the fixed degrees of freedom.
       * @return Number of values of a free row and a fixed column.
       */
      template <class Indices, class Fixed>
      static size_t countEliminated(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Fixed& fixed)
      {
        size_t res = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::countEliminated",
          Range(execution, 0, count), EliminatedCount<Indices, Fixed>{rows, cols, fixed},
          res);
        return res;
      }

      /**
       * @brief Clears entry e when its row is fixed and lists its value when
       * its column alone is fixed, at the scanned position.
       */
      template <class Indices, class Fixed, class Values, class Scalars>
      struct Elimination
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Fixed fixed; ///< Flags of the fixed degrees of freedom.
        Values values; ///< Fixed values, per degree of freedom.
        Scalars coo; ///< Values of the coordinate list.

        /// @brief Receives the rows of the moved values.
        ::Kokkos::View<Index*, MemorySpace> resRows;
        ::Kokkos::View<Real*, MemorySpace> resValues; ///< Receives the moved values.
        size_t offset; ///< Position of the first reserved contribution.

        /**
         * @brief Eliminates the fixed degrees of freedom from the entry.
         * @param e Position of the entry in the list.
         * @param update Number of values moved by the previous entries,
         * incremented.
         * @param final Whether the pass writes the results.
         */
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
       * @param execution Execution space instance of the kernel.
       * @param count Number of entries in the list.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param fixed Flags of the fixed degrees of freedom.
       * @param values Fixed values, per degree of freedom.
       * @param coo Values of the coordinate list, cleared on the fixed rows.
       * @param res Receives the moved values, in coordinate order.
       */
      template <class Indices, class Fixed, class Values, class Scalars>
      static void eliminate(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Fixed& fixed,
        const Values& values, const Scalars& coo, const Contributions& res)
      {
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::eliminate",
          Range(execution, 0, count),
          Elimination<Indices, Fixed, Values, Scalars>{
            rows, cols, fixed, values, coo, res.rows, res.values, res.offset});
      }

      /**
       * @brief Lists the entries of a free row and a fixed column, or those of
       * a fixed row, at their scanned positions.
       */
      template <class Indices, class Fixed>
      struct EntrySelection
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Fixed fixed; ///< Flags of the fixed degrees of freedom.
        bool
          moved; ///< Whether the entries of a free row and a fixed column are listed.
        ::Kokkos::View<Index*, MemorySpace>
          res; ///< Receives the entries, empty to count them.

        /**
         * @brief Lists the entry when it is selected.
         * @param e Position of the entry in the list.
         * @param update Number of selected entries before @p e, incremented.
         * @param final Whether the pass writes the results.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, size_t& update, const bool final) const
        {
          const bool free = !fixed(rows(e));
          if (moved ? !(free && fixed(cols(e))) : free)
            return;
          if (final && res.extent(0) > 0)
            res(update) = static_cast<Index>(e);
          ++update;
        }
      };

      /**
       * @brief Lists by increasing position the entries of a free row and a
       * fixed column, or those of a fixed row.
       * @param execution Execution space instance of the kernels.
       * @param count Number of entries in the list.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param fixed Flags of the fixed degrees of freedom.
       * @param moved Whether the entries whose value moves to the right-hand
       * side are listed, rather than those of a fixed row.
       * @return Positions of the selected entries.
       */
      template <class Indices, class Fixed>
      static ::Kokkos::View<Index*, MemorySpace> selectEntries(
        const ExecutionSpace& execution, size_t count, const Indices& rows,
        const Indices& cols, const Fixed& fixed, bool moved)
      {
        using Selection = EntrySelection<Indices, Fixed>;
        size_t selected = 0;
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::countSelected",
          Range(execution, 0, count), Selection{rows, cols, fixed, moved, {}}, selected);
        ::Kokkos::View<Index*, MemorySpace> res(
          ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
            std::string("Rodin eliminated entries")),
          selected);
        if (selected > 0)
        {
          ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::select",
            Range(execution, 0, count), Selection{rows, cols, fixed, moved, res});
        }
        return res;
      }

      /**
       * @brief Moves the value of a listed entry of a free row and a fixed
       * column to the right-hand side.
       */
      template <class Indices, class Values, class Scalars>
      struct MovedElimination
      {
        ::Kokkos::View<Index*, MemorySpace>
          entries; ///< Entries, by increasing position.
        Indices rows; ///< Row of each listed entry.
        Indices cols; ///< Column of each listed entry.
        Values values; ///< Fixed values, per degree of freedom.
        Scalars coo; ///< Values of the coordinate list.

        /// @brief Receives the rows of the moved values.
        ::Kokkos::View<Index*, MemorySpace> resRows;
        ::Kokkos::View<Real*, MemorySpace> resValues; ///< Receives the moved values.
        size_t offset; ///< Position of the first reserved contribution.

        /**
         * @brief Moves the value of an entry to the right-hand side.
         * @param k Rank of the entry among the listed ones.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t k) const
        {
          const Index e = entries(k);
          resRows(offset + k) = rows(k);
          resValues(offset + k) = -(coo(e) * values(cols(k)));
          coo(e) = 0;
        }
      };

      /// @brief Copies the coordinates of a listed entry.
      template <class Indices, class Results>
      struct EntryCoordinates
      {
        ::Kokkos::View<Index*, MemorySpace> entries; ///< Listed entries.
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Results resRows; ///< Receives the row of each listed entry.
        Results resCols; ///< Receives the column of each listed entry.

        /**
         * @brief Copies the coordinates of an entry.
         * @param k Rank of the entry among the listed ones.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t k) const
        {
          const Index e = entries(k);
          resRows(k) = static_cast<typename Results::value_type>(rows(e));
          resCols(k) = static_cast<typename Results::value_type>(cols(e));
        }
      };

      /// @brief Clears a listed entry of a fixed row.
      template <class Scalars>
      struct ClearedElimination
      {
        ::Kokkos::View<Index*, MemorySpace> entries; ///< Entries of a fixed row.
        Scalars coo; ///< Values of the coordinate list.

        /**
         * @brief Clears an entry.
         * @param k Rank of the entry among the listed ones.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t k) const
        {
          coo(entries(k)) = 0;
        }
      };

      /**
       * @brief Number of contributions of an entry of a linear integrator.
       * @param e Position of the entry.
       * @param indices Indices of the entries in the space of the integrator.
       * @param values Values of the entries.
       * @param row Offset of the space of the integrator.
       * @param hasFixed Whether the system has fixed degrees of freedom.
       * @param fixed Flags of the fixed degrees of freedom.
       * @param x Expansions on the masters, used when @p expanded is set.
       * @param expanded Whether an entry goes to the masters of its row, with
       * their coefficients, rather than to its row when it is free.
       * @return Zero for a zero value, the number of masters of the row when
       * @p expanded is set, and otherwise one when the row is free.
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
        Indices indices; ///< Indices of the entries in the space of the integrator.
        Values values; ///< Values of the entries.
        Index row; ///< Offset of the space of the integrator.
        bool hasFixed; ///< Whether the system has fixed degrees of freedom.
        Fixed fixed; ///< Flags of the fixed degrees of freedom.
        Expansions x; ///< Expansions on the masters.
        bool expanded; ///< Whether the entries go to the masters of their rows.

        /**
         * @brief Counts the contributions of the entry.
         * @param e Position of the entry.
         * @param sum Number of contributions, incremented.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += countEntry(e, indices, values, row, hasFixed, fixed, x, expanded);
        }
      };

      /**
       * @brief Counts the contributions of the entries of a linear
       * integrator to the right-hand side.
       * @param execution Execution space instance of the kernel.
       * @param begin Position of the first entry. @param count Entry count.
       * @param indices Indices of the entries in the space of the integrator.
       * @param values Values of the entries.
       * @param row Offset of the space of the integrator.
       * @param hasFixed Whether the system has fixed degrees of freedom.
       * @param fixed Flags of the fixed degrees of freedom.
       * @param x Expansions on the masters, used when @p expanded is set.
       * @param expanded Whether an entry goes to the masters of its row.
       * @return Number of contributions, as counted by countEntry.
       */
      template <class Indices, class Values, class Fixed>
      static size_t countVector(const ExecutionSpace& execution, size_t begin,
        size_t count, const Indices& indices, const Values& values, Index row,
        bool hasFixed, const Fixed& fixed, const Expansions& x, bool expanded)
      {
        size_t res = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::countVector",
          Range(execution, begin, begin + count),
          VectorCount<Indices, Values, Fixed>{
            indices, values, row, hasFixed, fixed, x, expanded},
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
        Indices indices; ///< Indices of the entries in the space of the integrator.
        Values values; ///< Values of the entries.
        Index row; ///< Offset of the space of the integrator.
        bool hasFixed; ///< Whether the system has fixed degrees of freedom.
        Fixed fixed; ///< Flags of the fixed degrees of freedom.
        Expansions x; ///< Expansions on the masters.
        bool expanded; ///< Whether the entries go to the masters of their rows.

        /// @brief Receives the rows of the contributions.
        ::Kokkos::View<Index*, MemorySpace> resRows;
        ::Kokkos::View<Real*, MemorySpace> resValues; ///< Receives the contributions.
        size_t offset; ///< Position of the first reserved contribution.

        /**
         * @brief Lists the contributions of the entry.
         * @param e Position of the entry.
         * @param update Number of contributions of the previous entries,
         * incremented.
         * @param final Whether the pass writes the results.
         */
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
       * @param execution Execution space instance of the kernel.
       * @param begin Position of the first entry. @param count Entry count.
       * @param indices Indices of the entries in the space of the integrator.
       * @param values Values of the entries.
       * @param row Offset of the space of the integrator.
       * @param hasFixed Whether the system has fixed degrees of freedom.
       * @param fixed Flags of the fixed degrees of freedom.
       * @param x Expansions on the masters, used when @p expanded is set.
       * @param expanded Whether an entry goes to the masters of its row.
       * @param res Receives the contributions, counted by countVector.
       *
       * The contributions follow the order of the entries, and of the masters
       * of an entry.
       */
      template <class Indices, class Values, class Fixed>
      static void listVector(const ExecutionSpace& execution, size_t begin, size_t count,
        const Indices& indices, const Values& values, Index row, bool hasFixed,
        const Fixed& fixed, const Expansions& x, bool expanded, const Contributions& res)
      {
        const auto policy = Range(execution, begin, begin + count);
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::vector", policy,
          VectorList<Indices, Values, Fixed>{indices, values, row, hasFixed, fixed, x,
            expanded, res.rows, res.values, res.offset});
      }

      /// @brief Writes the fixed value of row i when it is fixed.
      template <class Fixed, class Values, class Owned, class RHS>
      struct FixedImposition
      {
        Fixed fixed; ///< Flags of the rows whose value is imposed.
        Values values; ///< Imposed values, per degree of freedom.

        /**
         * @brief Flags of the degrees of freedom the process owns, empty when the system
         * is not distributed.
         */
        Owned owned;
        RHS rhs; ///< Right-hand side.
        bool distributed; ///< Whether the system is distributed.

        /**
         * @brief Writes the fixed value of the row.
         * @param i Row.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          if (fixed(i))
            rhs(i) = !distributed || owned(i) ? values(i) : Real(0);
        }
      };

      /**
       * @brief Writes the fixed values on the fixed rows.
       * @param execution Execution space instance of the kernel.
       * @param count Number of rows.
       * @param fixed Flags of the rows whose value is imposed.
       * @param values Imposed values, per degree of freedom.
       * @param owned Flags of the degrees of freedom the process owns, empty
       * when the system is not distributed.
       * @param rhs Right-hand side.
       *
       * On a distributed system, the other processes write zero.
       */
      template <class Fixed, class Values, class Owned, class RHS>
      static void imposeFixed(const ExecutionSpace& execution, size_t count,
        const Fixed& fixed, const Values& values, const Owned& owned, const RHS& rhs)
      {
        const bool distributed = owned.extent(0) > 0;
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::fixed",
          Range(execution, 0, count),
          FixedImposition<Fixed, Values, Owned, RHS>{
            fixed, values, owned, rhs, distributed});
      }

      /// @brief Writes the row of entry e of a linear integrator.
      template <class Indices>
      struct EntryRows
      {
        Indices indices; ///< Indices of the entries in the space of the integrator.
        Index row; ///< Offset of the test space.
        Indices res; ///< Receives the rows.

        /**
         * @brief Writes the row of the entry.
         * @param e Position of the entry.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          res(e) = row + indices(e);
        }
      };

      /**
       * @brief Writes the rows of the entries of a linear integrator.
       * @param execution Execution space instance of the kernel.
       * @param indices Entry indices. @param begin Position of the first
       * entry. @param end Position past the last one.
       * @param row Offset of the test space. @param res Receives the rows.
       */
      template <class Indices>
      static void shiftEntries(const ExecutionSpace& execution, const Indices& indices,
        size_t begin, size_t end, Index row, const Indices& res)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::rows",
          Range(execution, begin, end), EntryRows<Indices>{indices, row, res});
      }

      /// @brief Numbers entry i of a dense vector.
      template <class Indices>
      struct Sequence
      {
        Indices res; ///< Receives the numbering.

        /**
         * @brief Writes the index of the entry into it.
         * @param i Index of the entry.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          res(i) = i;
        }
      };

      /**
       * @brief Numbers the entries of a dense vector.
       * @param execution Execution space instance of the kernel.
       * @param res Receives the index of each entry.
       */
      template <class Indices>
      static void sequence(const ExecutionSpace& execution, const Indices& res)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::sequence",
          Range(execution, 0, res.extent(0)), Sequence<Indices>{res});
      }

      /**
       * @brief Counts the contributions of entry e expanded on the masters,
       * and the affine values it lifts.
       */
      template <class Indices, class Scalars>
      struct ExpandedCount
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Scalars values; ///< Values of the coordinate list.
        Expansions x; ///< Expansions on the masters.

        /**
         * @brief Counts the contributions and lifted values of the entry.
         * @param e Position of the entry in the list.
         * @param sum Numbers of contributions and of lifted values, incremented.
         */
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
       * @param execution Execution space instance of the kernel.
       * @param count Number of entries in the list.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param values Values of the coordinate list.
       * @param x Expansions on the masters.
       * @return Numbers of expanded contributions and of lifted values.
       */
      template <class Indices, class Scalars>
      static Sizes countExpanded(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Scalars& values,
        const Expansions& x)
      {
        Sizes res;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::countExpanded",
          Range(execution, 0, count),
          ExpandedCount<Indices, Scalars>{rows, cols, values, x}, res);
        return res;
      }

      /**
       * @brief Expands entry e on the masters and lists the affine values it
       * lifts, at the scanned positions.
       */
      template <class Indices, class Results, class Scalars>
      struct Expansion
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Scalars values; ///< Values of the coordinate list.
        Expansions x; ///< Expansions on the masters.
        Results resRows; ///< Receives the rows of the expanded contributions.
        Results resCols; ///< Receives the columns of the expanded contributions.
        Scalars resValues; ///< Receives the expanded contributions.

        /// @brief Receives the rows of the lifted values.
        ::Kokkos::View<Index*, MemorySpace> liftRows;
        ::Kokkos::View<Real*, MemorySpace> liftValues; ///< Receives the lifted values.
        size_t offset; ///< Position of the first reserved lifted value.

        /**
         * @brief Expands the entry and lists the values it lifts.
         * @param e Position of the entry in the list.
         * @param update Numbers of contributions and of lifted values of the
         * previous entries, incremented.
         * @param final Whether the pass writes the results.
         */
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
       * @param execution Execution space instance of the kernel.
       * @param count Number of entries in the list.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param values Values of the coordinate list.
       * @param x Expansions on the masters.
       * @param resRows Receives the rows of the expanded contributions.
       * @param resCols Receives the columns of the expanded contributions.
       * @param resValues Receives the expanded contributions.
       * @param lifts Receives the lifted affine values, in coordinate order.
       *
       * The nonzero value @f$ a @f$ of row @f$ i @f$ and column @f$ j @f$
       * becomes the values @f$ (c_r a) c_s @f$ of the masters @f$ r @f$ of
       * @f$ i @f$ and @f$ s @f$ of @f$ j @f$, and lifts @f$ -(c_r a) d_j @f$
       * to row @f$ r @f$ when @f$ j @f$ has a nonzero affine value
       * @f$ d_j @f$, as the host assemblies do.
       */
      template <class Indices, class Results, class Scalars>
      static void expand(const ExecutionSpace& execution, size_t count,
        const Indices& rows, const Indices& cols, const Scalars& values,
        const Expansions& x, const Results& resRows, const Results& resCols,
        const Scalars& resValues, const Contributions& lifts)
      {
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::expand",
          Range(execution, 0, count),
          Expansion<Indices, Results, Scalars>{rows, cols, values, x, resRows, resCols,
            resValues, lifts.rows, lifts.values, lifts.offset});
      }

      /// @brief Counts position i when two index lists differ there.
      template <class Indices, class Previous>
      struct DifferenceCount
      {
        Indices indices; ///< Index list.
        Previous previous; ///< Index list it is compared with.

        /**
         * @brief Counts the position when the lists differ there.
         * @param i Position in the lists.
         * @param sum Number of differing positions, incremented.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i, size_t& sum) const
        {
          sum += indices(i) != previous(i);
        }
      };

      /**
       * @brief Counts the positions at which two index lists differ.
       * @param execution Execution space instance of the kernel.
       * @param indices Index list.
       * @param previous Index list it is compared with, at least as long.
       * @return Number of differing positions.
       */
      template <class Indices, class Previous>
      static size_t countDifferences(
        const ExecutionSpace& execution, const Indices& indices, const Previous& previous)
      {
        size_t differences = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::differences",
          Range(execution, 0, indices.extent(0)),
          DifferenceCount<Indices, Previous>{indices, previous}, differences);
        return differences;
      }

      /// @brief Copies coordinate i into another index type.
      template <class From, class To>
      struct CoordinateCopy
      {
        From from; ///< Coordinates copied.
        To to; ///< Receives the coordinates.
        size_t offset; ///< Position of the first coordinate received.

        /**
         * @brief Copies the coordinate.
         * @param i Position of the coordinate copied.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          to(offset + i) = static_cast<typename To::value_type>(from(i));
        }
      };

      /**
       * @brief Copies coordinates into another index type.
       * @param execution Execution space instance of the kernel.
       * @param from Coordinates copied.
       * @param to Receives the coordinates.
       * @param offset Position of the first coordinate received.
       */
      template <class From, class To>
      static void copyCoordinates(const ExecutionSpace& execution, const From& from,
        const To& to, size_t offset = 0)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::copyCoordinates",
          Range(execution, 0, from.extent(0)),
          CoordinateCopy<From, To>{from, to, offset});
      }

      /// @brief Counts coordinate e when it differs between two patterns.
      template <class Indices, class Previous>
      struct ChangeCount
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Previous previousRows; ///< Row coordinates of the previous pattern.
        Previous previousCols; ///< Column coordinates of the previous pattern.

        /**
         * @brief Counts the coordinate when it differs between the patterns.
         * @param e Position of the coordinate.
         * @param sum Number of differing coordinates, incremented.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += static_cast<Index>(rows(e)) != static_cast<Index>(previousRows(e)) ||
            static_cast<Index>(cols(e)) != static_cast<Index>(previousCols(e));
        }
      };

      /**
       * @brief Counts the coordinates that differ between two patterns.
       * @param execution Execution space instance of the kernel.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param previousRows Row coordinates of the previous pattern.
       * @param previousCols Column coordinates of the previous pattern.
       * @param begin Position of the first coordinate compared.
       * @param end Position past the last coordinate compared, the end of
       * the list by default.
       * @return Number of differing coordinates.
       */
      template <class Indices, class Previous>
      static size_t countChanges(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, const Previous& previousRows, const Previous& previousCols,
        size_t begin = 0, size_t end = std::numeric_limits<size_t>::max())
      {
        size_t changes = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::pattern",
          Range(execution, begin, std::min<size_t>(end, rows.extent(0))),
          ChangeCount<Indices, Previous>{rows, cols, previousRows, previousCols},
          changes);
        return changes;
      }

      /**
       * @brief Sorts the coordinates by row then column.
       * @param execution Execution space instance of the kernel.
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

      /// @brief Coordinate followed by the original contribution rank.
      struct PatternKey
      {
        std::uint64_t coordinate; ///< Row times column count plus column.
        size_t rank; ///< Rank in the concatenated integrator listings.
      };

      /// @brief Total order, including a deterministic order of duplicates.
      struct PatternLess
      {
        /**
         * @brief Orders coordinates, breaking ties by contribution rank.
         * @param lhs Left key.
         * @param rhs Right key.
         * @returns Whether the left key precedes the right key.
         */
        KOKKOS_INLINE_FUNCTION bool operator()(
          const PatternKey& lhs, const PatternKey& rhs) const
        {
          return lhs.coordinate < rhs.coordinate ||
            (lhs.coordinate == rhs.coordinate && lhs.rank < rhs.rank);
        }
      };

      /**
       * @brief Largest number of coordinates that sortPattern() sorts at
       * once.
       *
       * The key array is limited to the same one-GiB budget as the former
       * eight-byte keys. Each key now also stores the contribution rank.
       * Row ranges are halved until they fit this budget; an indivisible
       * single row may exceed it. Sort scratch is proportional to the range.
       */
      static constexpr size_t SortBlockSize = (size_t(1) << 30) / sizeof(PatternKey);

      /// @brief Finds the smallest and the largest row coordinates.
      template <class Indices>
      struct RowBounds
      {
        Indices rows; ///< Row coordinates.

        /**
         * @brief Folds the row of a coordinate into the bounds.
         * @param e Index of the coordinate.
         * @param smallest Smallest row so far.
         * @param largest Largest row so far.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, Index& smallest, Index& largest) const
        {
          const Index row = rows(e);
          smallest = row < smallest ? row : smallest;
          largest = row > largest ? row : largest;
        }
      };

      /// @brief Counts the coordinates whose row lies in a range.
      template <class Indices>
      struct RowCount
      {
        Indices rows; ///< Row coordinates.
        Index low; ///< First row of the range.
        Index high; ///< Row past the range.

        /**
         * @brief Counts the coordinate when its row lies in the range.
         * @param e Index of the coordinate.
         * @param count Number of coordinates in the range, incremented.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& count) const
        {
          const Index row = rows(e);
          if (row >= low && row < high)
            ++count;
        }
      };

      /**
       * @brief Lists, in rank order, the sort key and the index of the
       * coordinates of a segment whose row lies in a range.
       */
      template <class Indices, class Order>
      struct PatternKeys
      {
        /// @brief Type of the indices in the permutation.
        using Position = typename Order::value_type;

        size_t begin; ///< Position of the first entry of the segment.
        size_t size; ///< Number of entries of the segment.
        size_t spacing; ///< Spacing of the entries of the segment.
        Index low; ///< First row of the range.
        Index high; ///< Row past the range.
        size_t base; ///< Number of coordinates of the range listed before the segment.

        /// @brief Receives the keys, the row times the column count plus the column.
        ::Kokkos::View<PatternKey*, MemorySpace> keys;
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        size_t columns; ///< Column count of the system.
        Order order; ///< Receives the indices of the coordinates of the range.

        /**
         * @brief Lists the coordinate of rank r when its row lies in the
         * range.
         * @param r Rank of the coordinate in the segment.
         * @param listed Number of coordinates of the range listed before it,
         * incremented.
         * @param final Whether the pass writes the results.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t r, size_t& listed, const bool final) const
        {
          const size_t e =
            begin + ::Rodin::Kokkos::Variational::getEntryPosition(r, size, spacing);
          const Index row = rows(e);
          if (row < low || row >= high)
            return;
          if (final)
          {
            keys(base + listed) = {std::uint64_t(row) * columns + cols(e), begin + r};
            order(base + listed) = static_cast<Position>(e);
          }
          ++listed;
        }
      };

      /**
       * @brief Sorts the coordinates into the given indices.
       * @param execution Execution space instance of the kernel.
       * @param columns Column count of the system.
       * @param segments Segments of the coordinates.
       * @param offsets Position of the first entry of every segment.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param order Receives the indices of the sorted coordinates.
       * @param blockSize Target maximum number of coordinates sorted at once.
       *
       * The rows are split in ranges of at most SortBlockSize coordinates,
       * by halving, and the ranges are sorted one after the other in row
       * order. Contribution rank breaks coordinate ties, so every bucket
       * size and backend produces the same permutation.
       */
      template <class Indices, class Order>
      static void sortPattern(const ExecutionSpace& execution, size_t columns,
        const std::vector<COOSegment>& segments, const std::vector<size_t>& offsets,
        const Indices& rows, const Indices& cols, const Order& order,
        size_t blockSize = SortBlockSize)
      {
        const size_t count = rows.extent(0);
        if (count == 0)
          return;
        Index smallest = 0;
        Index largest = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::rowBounds",
          Range(execution, 0, count), RowBounds<Indices>{rows},
          ::Kokkos::Min<Index>(smallest), ::Kokkos::Max<Index>(largest));
        struct Rows
        {
          Index low;
          Index high;
          size_t count;
        };
        std::vector<Rows> pending{{smallest, largest + 1, count}};
        ::Kokkos::View<PatternKey*, MemorySpace> keys;
        size_t sorted = 0;
        while (!pending.empty())
        {
          const Rows range = pending.back();
          pending.pop_back();
          if (range.count == 0)
            continue;
          if (range.count > blockSize && range.high - range.low > 1)
          {
            const Index middle = range.low + (range.high - range.low) / 2;
            size_t below = 0;
            ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::rowCount",
              Range(execution, 0, count), RowCount<Indices>{rows, range.low, middle},
              below);
            pending.push_back({middle, range.high, range.count - below});
            pending.push_back({range.low, middle, below});
            continue;
          }
          if (keys.extent(0) < range.count)
          {
            keys = {};
            keys = ::Kokkos::View<PatternKey*, MemorySpace>(
              ::Kokkos::view_alloc(
                execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO keys")),
              std::max(range.count, std::min(count, blockSize)));
          }
          const auto block =
            ::Kokkos::subview(order, std::make_pair(sorted, sorted + range.count));
          size_t listed = 0;
          for (size_t s = 0; s <= segments.size(); ++s)
          {
            const size_t begin = offsets[s];
            const size_t size = (s < segments.size() ? offsets[s + 1] : count) - begin;
            const size_t spacing = s < segments.size() ? segments[s].spacing : 1;
            size_t found = 0;
            ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::keys",
              Range(execution, 0, size),
              PatternKeys<Indices, decltype(block)>{begin, size, spacing, range.low,
                range.high, listed, keys, rows, cols, columns, block},
              found);
            listed += found;
          }
          ::Kokkos::Experimental::sort_by_key(execution,
            ::Kokkos::subview(keys, std::make_pair(size_t(0), range.count)), block,
            PatternLess{});
          sorted += range.count;
        }
      }

      /// @brief Whether the sorted coordinate i differs from the previous one.
      template <class Indices, class Order>
      struct DistinctHead
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Order order; ///< Indices of the sorted coordinates.

        /**
         * @brief Tests the sorted coordinate.
         * @param i Position in the sorted coordinates.
         * @return Whether the coordinate is the first or differs from the
         * previous one.
         */
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
        /// @brief Test of the starts of the distinct coordinates.
        DistinctHead<Indices, Order> head;

        /**
         * @brief Counts the sorted coordinate when it starts a distinct one.
         * @param i Position in the sorted coordinates.
         * @param sum Number of distinct coordinates, incremented.
         */
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
        /// @brief Test of the starts of the distinct coordinates.
        DistinctHead<Indices, Order> head;

        /// @brief Receives the start of each distinct coordinate.
        ::Kokkos::View<Index*, MemorySpace> starts;

        /**
         * @brief Writes the position of the sorted coordinate when it starts a
         * distinct one.
         * @param i Position in the sorted coordinates.
         * @param update Number of distinct coordinates before it, incremented.
         * @param final Whether the pass writes the results.
         */
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
       * @param execution Execution space instance of the kernel.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param order Indices of the sorted coordinates.
       * @returns Starts of the distinct coordinates, followed by the list size.
       */
      template <class Indices, class Order>
      static ::Kokkos::View<Index*, MemorySpace> distinctStarts(
        const ExecutionSpace& execution, const Indices& rows, const Indices& cols,
        const Order& order)
      {
        const size_t count = order.extent(0);
        const DistinctHead<Indices, Order> head{rows, cols, order};
        size_t distinct = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::PETSc::distinct",
          Range(execution, 0, count), DistinctCount<Indices, Order>{head}, distinct);
        ::Kokkos::View<Index*, MemorySpace> starts(
          ::Kokkos::view_alloc(
            execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO starts")),
          distinct + 1);
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::PETSc::starts",
          Range(execution, 0, count), DistinctStartList<Indices, Order>{head, starts});
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(starts, distinct), count);
        return starts;
      }

      /// @brief Lists the distinct coordinate u in the PETSc index type.
      template <class Indices, class Order, class Coordinates>
      struct DistinctList
      {
        Indices rows; ///< Row coordinates.
        Indices cols; ///< Column coordinates.
        Order order; ///< Indices of the sorted coordinates.

        /// @brief Start of each distinct coordinate in the sorted coordinates.
        ::Kokkos::View<Index*, MemorySpace> starts;
        Coordinates distinctRows; ///< Receives the rows of the distinct coordinates.
        Coordinates distinctCols; ///< Receives the columns of the distinct coordinates.

        /**
         * @brief Lists the distinct coordinate.
         * @param u Index of the distinct coordinate.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          const Index e = order(starts(u));
          distinctRows(u) = static_cast<PetscInt>(rows(e));
          distinctCols(u) = static_cast<PetscInt>(cols(e));
        }
      };

      /**
       * @brief Lists the distinct coordinates in the PETSc index type.
       * @param execution Execution space instance of the kernel.
       * @param rows Row coordinates.
       * @param cols Column coordinates.
       * @param order Indices of the sorted coordinates.
       * @param starts Start of each distinct coordinate in the sorted
       * coordinates.
       * @param distinctRows Receives the rows of the distinct coordinates.
       * @param distinctCols Receives the columns of the distinct coordinates.
       */
      template <class Indices, class Order, class Coordinates>
      static void listDistinct(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, const Order& order,
        const ::Kokkos::View<Index*, MemorySpace>& starts,
        const Coordinates& distinctRows, const Coordinates& distinctCols)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::distinctCoordinates",
          Range(execution, 0, distinctRows.extent(0)),
          DistinctList<Indices, Order, Coordinates>{
            rows, cols, order, starts, distinctRows, distinctCols});
      }

      /// @brief Writes the distinct coordinate u at its entries of the list.
      template <class Indices, class Order, class Coordinates>
      struct Restoration
      {
        Order order; ///< Indices of the sorted coordinates.

        /// @brief Start of each distinct coordinate in the sorted coordinates.
        ::Kokkos::View<Index*, MemorySpace> starts;
        Coordinates distinctRows; ///< Rows of the distinct coordinates.
        Coordinates distinctCols; ///< Columns of the distinct coordinates.
        Indices rows; ///< Receives the row coordinates.
        Indices cols; ///< Receives the column coordinates.

        /**
         * @brief Writes the distinct coordinate at its entries.
         * @param u Index of the distinct coordinate.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          const auto row = static_cast<typename Indices::value_type>(distinctRows(u));
          const auto col = static_cast<typename Indices::value_type>(distinctCols(u));
          for (Index i = starts(u); i < starts(u + 1); ++i)
          {
            const Index e = order(i);
            rows(e) = row;
            cols(e) = col;
          }
        }
      };

      /**
       * @brief Lists the coordinates again from the distinct ones, the
       * inverse of listDistinct().
       * @param execution Execution space instance of the kernel.
       * @param order Indices of the sorted coordinates.
       * @param starts Start of each distinct coordinate in the sorted
       * coordinates.
       * @param distinctRows Rows of the distinct coordinates.
       * @param distinctCols Columns of the distinct coordinates.
       * @param rows Receives the row coordinates.
       * @param cols Receives the column coordinates.
       */
      template <class Indices, class Order, class Coordinates>
      static void restore(const ExecutionSpace& execution, const Order& order,
        const ::Kokkos::View<Index*, MemorySpace>& starts,
        const Coordinates& distinctRows, const Coordinates& distinctCols,
        const Indices& rows, const Indices& cols)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::restoreCoordinates",
          Range(execution, 0, distinctRows.extent(0)),
          Restoration<Indices, Order, Coordinates>{
            order, starts, distinctRows, distinctCols, rows, cols});
      }

      /// @brief Sums the values of the distinct coordinate u.
      template <class Order, class Indices, class Scalars>
      struct Merge
      {
        Order order; ///< Indices of the sorted coordinates.

        /**
         * @brief Start of each distinct coordinate in the sorted coordinates, followed
         * by the list size.
         */
        Indices starts;
        Scalars values; ///< Values of the coordinate list.
        Scalars merged; ///< Receives the sum of each distinct coordinate.

        /**
         * @brief Sums the values of the distinct coordinate.
         * @param u Index of the distinct coordinate.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          PetscScalar sum = 0;
          for (Index i = starts(u); i < starts(u + 1); ++i)
            sum += values(order(i));
          merged(u) = sum;
        }
      };

      /**
       * @brief Sums the values of each distinct coordinate.
       * @param execution Execution space instance of the kernel.
       * @param order Indices of the sorted coordinates.
       * @param starts Start of each distinct coordinate in the sorted
       * coordinates, followed by the list size.
       * @param values Values of the coordinate list.
       * @param merged Receives the sum of each distinct coordinate.
       */
      template <class Order, class Indices, class Scalars>
      static void merge(const ExecutionSpace& execution, const Order& order,
        const Indices& starts, const Scalars& values, const Scalars& merged)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::merge",
          Range(execution, 0, merged.extent(0)),
          Merge<Order, Indices, Scalars>{order, starts, values, merged});
      }

    private:
      struct Matrix
      {
        BilinearType* integrator;
        Index row;
        Index col;
        Real sign;
      };

      struct Vector
      {
        LinearType* integrator;
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
       * @return Coordinates and values of the entries.
       */
      Entries listOperators() const
      {
        Entries res;
        std::unordered_map<PetscInt, Index> locals;
        if (m_distribution && !m_operators.empty())
          locals = m_distribution->getLocals();
        const auto local =
          [&](PetscInt global) {
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
          const PetscInt rowOffset =
            m_distribution ? m_distribution->getOffset(o.row) : 0;
          const PetscInt colOffset =
            m_distribution ? m_distribution->getOffset(o.col) : 0;
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

      /**
       * @brief Unmanaged host view of a vector.
       * @param data Vector viewed, alive while the view is used.
       * @return View of the values of the vector.
       */
      template <class T>
      static auto host(const std::vector<T>& data)
      {
        return ::Kokkos::View<const T*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(data.data(), data.size());
      }

      /**
       * @brief Reallocates a view, without initializing it, when its size
       * changes.
       * @param view View resized, kept when it already has the size.
       * @param size Number of entries.
       * @param label Label of the allocation.
       */
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

      /**
       * @brief Copies a host vector to a view, resized to its size.
       * @param view Destination view.
       * @param host Values copied.
       * @param label Label of the allocation, when the view is resized.
       */
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
       * @brief Keeps the fixed degrees of freedom of the assembly, whose
       * revision increases when they change.
       * @param cache Buffers and pattern of the previous assemblies.
       * @param mask Flags of the fixed degrees of freedom.
       */
      static void keepFixed(CacheType& cache, const std::vector<unsigned char>& mask)
      {
        if (cache.fixedMask == mask)
          return;
        cache.fixedMask = mask;
        ++cache.fixedRevision;
      }

      /**
       * @brief Releases the coordinates of the sorted pattern, keeping the
       * distinct ones.
       *
       * The coordinates take half of the memory of the list, and the
       * assemblies that keep the pattern do not read them.
       * restoreCoordinates() lists them again for the ones that do.
       * @param cache Buffers and pattern of the previous assemblies.
       */
      void releaseCoordinates(CacheType& cache) const
      {
        if (cache.released || !cache.indexed)
          return;
        const size_t distinct = cache.starts.extent(0) - 1;
        cache.distinctRows = decltype(cache.distinctRows)(
          ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
            std::string("Rodin distinct COO rows")),
          distinct);
        cache.distinctCols = decltype(cache.distinctCols)(
          ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
            std::string("Rodin distinct COO cols")),
          distinct);
        cache.order.visit([&](const auto& order) {
          listDistinct(m_execution, cache.rows, cache.cols, order, cache.starts,
            cache.distinctRows, cache.distinctCols);
        });
        cache.rows = {};
        cache.cols = {};
        cache.released = true;
      }

      /**
       * @brief Lists again the coordinates released by releaseCoordinates().
       * @param cache Buffers and pattern of the previous assemblies.
       */
      void restoreCoordinates(CacheType& cache) const
      {
        if (!cache.released)
          return;
        cache.rows = decltype(cache.rows)(
          ::Kokkos::view_alloc(
            m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO rows")),
          cache.listed);
        cache.cols = decltype(cache.cols)(
          ::Kokkos::view_alloc(
            m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO cols")),
          cache.listed);
        cache.order.visit([&](const auto& order) {
          restore(m_execution, order, cache.starts, cache.distinctRows,
            cache.distinctCols, cache.rows, cache.cols);
        });
        cache.distinctRows = {};
        cache.distinctCols = {};
        cache.released = false;
      }

      /**
       * @brief Integrates the local integrators into the coordinate list.
       * @param cache Buffers and pattern of the previous assemblies.
       * @param fixed Whether the fixed degrees of freedom are eliminated.
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

        // When only the restricted segments may have moved, they alone are
        // listed again and compared with the kept pattern.
        const bool kept = cache.indexed && m_operators.empty() &&
          cache.segments == segments && cache.listed == count;
        const bool index = !kept || restricted;
        const bool partial = kept && restricted;
        // The restricted segments are packed, with their values when they
        // are listed again, so that the listing takes their size only.
        std::vector<bool> movable(segments.size());
        std::vector<size_t> packed(segments.size() + 1, 0);
        for (size_t s = 0; s < segments.size(); ++s)
        {
          movable[s] = !segments[s].layout.attributes.empty() ||
            !segments[s].layout.trialAttributes.empty();
          packed[s + 1] = packed[s] + (movable[s] ? segments[s].count : 0);
        }
        const std::vector<bool> listed =
          partial ? movable : std::vector<bool>(segments.size(), index);
        resize(cache.values, count, "Rodin COO values");
        if (partial)
        {
          resize(cache.nextRows, packed.back(), "Rodin listed COO rows");
          resize(cache.nextCols, packed.back(), "Rodin listed COO cols");
          resize(cache.nextValues, packed.back(), "Rodin listed COO values");
        }
        else if (index)
        {
          resize(cache.nextRows, count, "Rodin listed COO rows");
          resize(cache.nextCols, count, "Rodin listed COO cols");
        }
        const ::Kokkos::View<Index*, MemorySpace> none;
        std::vector<bool> filled(segments.size(), false);
        ::Rodin::Kokkos::Internal::clearLaunches();
        for (size_t k = 0; k < m_matrices.size(); ++k)
        {
          const auto& m = m_matrices[k];
          const size_t s = segment[k];
          const bool first = !filled[s];
          const bool pack = partial && listed[s];
          m.integrator->fill(pack ? cache.nextValues : cache.values,
            pack ? packed[s] : offsets[s], m.sign, !first,
            first && listed[s] ? cache.nextRows : none,
            first && listed[s] ? cache.nextCols : none);
          ::Rodin::Kokkos::Internal::checkLaunches("A device bilinear integrator");
          filled[s] = true;
        }
        const auto range = std::make_pair(offsets.back(), tail);
        ::Kokkos::deep_copy(
          m_execution, ::Kokkos::subview(cache.values, range), host(entries.values));
        if (partial)
        {
          // The other segments, the diagonal and the kept pattern stay.
          bool changed = false;
          for (size_t s = 0; s < segments.size(); ++s)
          {
            if (!listed[s])
              continue;
            const auto moved = std::make_pair(offsets[s], offsets[s + 1]);
            const auto from = std::make_pair(packed[s], packed[s + 1]);
            ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(cache.values, moved),
              ::Kokkos::subview(cache.nextValues, from));
            shift(m_execution, cache.nextRows, cache.nextCols, packed[s],
              segments[s].count, segments[s].row, segments[s].col);
            const ::Kokkos::View<Index*, MemorySpace> nextRows =
              ::Kokkos::subview(cache.nextRows, from);
            const ::Kokkos::View<Index*, MemorySpace> nextCols =
              ::Kokkos::subview(cache.nextCols, from);
            const ::Kokkos::View<Index*, MemorySpace> rows =
              ::Kokkos::subview(cache.restrictedRows, from);
            const ::Kokkos::View<Index*, MemorySpace> cols =
              ::Kokkos::subview(cache.restrictedCols, from);
            if (countChanges(m_execution, nextRows, nextCols, rows, cols) == 0)
              continue;
            ::Kokkos::deep_copy(m_execution, rows, nextRows);
            ::Kokkos::deep_copy(m_execution, cols, nextCols);
            // The pattern is sorted again from all its coordinates.
            restoreCoordinates(cache);
            copyCoordinates(m_execution, nextRows, cache.rows, offsets[s]);
            copyCoordinates(m_execution, nextCols, cache.cols, offsets[s]);
            changed = true;
          }
          if (changed)
          {
            sortPattern(m_execution, m_cols, cache.segments, offsets, cache.rows,
              cache.cols, cache.order);
            cache.starts = cache.order.visit([&](const auto& order) {
              return distinctStarts(m_execution, cache.rows, cache.cols, order);
            });
            cache.matrix = -1;
            ++cache.pattern;
          }
        }
        else if (index)
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
          bool changed = !cache.indexed || cache.listed != count;
          if (!changed)
          {
            restoreCoordinates(cache);
            changed = countChanges(m_execution, cache.nextRows, cache.nextCols,
                        cache.rows, cache.cols) > 0;
          }
          resize(cache.restrictedRows, packed.back(), "Rodin restricted COO rows");
          resize(cache.restrictedCols, packed.back(), "Rodin restricted COO cols");
          for (size_t s = 0; s < segments.size(); ++s)
          {
            if (!movable[s])
              continue;
            const auto from = std::make_pair(offsets[s], offsets[s + 1]);
            const auto to = std::make_pair(packed[s], packed[s + 1]);
            ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(cache.restrictedRows, to),
              ::Kokkos::subview(cache.nextRows, from));
            ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(cache.restrictedCols, to),
              ::Kokkos::subview(cache.nextCols, from));
          }
          // The previous coordinates are released rather than kept for the
          // next listing, which is partial or rare, and the listed ones are
          // kept in the PETSc index type, one array after the other.
          assert(std::max(m_rows, m_cols) <=
            static_cast<size_t>(std::numeric_limits<PetscInt>::max()));
          cache.rows = {};
          cache.cols = {};
          cache.rows = decltype(cache.rows)(
            ::Kokkos::view_alloc(
              m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO rows")),
            count);
          copyCoordinates(m_execution, cache.nextRows, cache.rows);
          cache.nextRows = {};
          cache.cols = decltype(cache.cols)(
            ::Kokkos::view_alloc(
              m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO cols")),
            count);
          copyCoordinates(m_execution, cache.nextCols, cache.cols);
          cache.nextCols = {};
          cache.distinctRows = {};
          cache.distinctCols = {};
          cache.released = false;
          cache.listed = count;
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
            ++cache.pattern;
          }
        }
        fillDiagonal(
          m_execution, tail, diagonal, fixed, cache.fixed, m_owned, cache.values);
        if (fixed)
          eliminateListed(cache, tail);
      }

      /**
       * @brief Eliminates the fixed degrees of freedom from the listed
       * coordinates, as eliminate() does.
       *
       * The entries a fixed degree of freedom reaches, a small part of the
       * list, depend on the pattern and on the fixed degrees of freedom
       * alone. They are listed again when one of them changes, and only
       * they are visited at the other assemblies.
       * @param cache Buffers and pattern of the previous assemblies.
       * @param count Number of entries eliminated, at the start of the list.
       */
      void eliminateListed(CacheType& cache, size_t count) const
      {
        const std::pair key{cache.pattern, cache.fixedRevision};
        using Indices = decltype(cache.movedRows);
        if (cache.elimination != key || cache.eliminationCount != count)
        {
          restoreCoordinates(cache);
          cache.movedEntries =
            selectEntries(m_execution, count, cache.rows, cache.cols, cache.fixed, true);
          cache.clearedEntries =
            selectEntries(m_execution, count, cache.rows, cache.cols, cache.fixed, false);
          // The moved entries keep their coordinates, read at every assembly.
          resize(cache.movedRows, cache.movedEntries.extent(0), "Rodin eliminated rows");
          resize(cache.movedCols, cache.movedEntries.extent(0), "Rodin eliminated cols");
          ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::eliminatedCoordinates",
            Range(m_execution, 0, cache.movedEntries.extent(0)),
            EntryCoordinates<decltype(cache.rows), Indices>{cache.movedEntries,
              cache.rows, cache.cols, cache.movedRows, cache.movedCols});
          cache.elimination = key;
          cache.eliminationCount = count;
        }
        const size_t moved = cache.movedEntries.extent(0);
        const Contributions res = reserve(cache, moved);
        using Values = decltype(cache.fixedValues);
        using Scalars = decltype(cache.values);
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::eliminateMoved",
          Range(m_execution, 0, moved),
          MovedElimination<Indices, Values, Scalars>{cache.movedEntries, cache.movedRows,
            cache.movedCols, cache.fixedValues, cache.values, res.rows, res.values,
            res.offset});
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::PETSc::eliminateCleared",
          Range(m_execution, 0, cache.clearedEntries.extent(0)),
          ClearedElimination<Scalars>{cache.clearedEntries, cache.values});
      }

      /**
       * @brief Expands the coordinate list on the masters, appends the rows
       * of the identified degrees of freedom and every diagonal entry, then
       * eliminates the fixed degrees of freedom.
       * @param cache Buffers and pattern of the previous assemblies, holding
       * the listed coordinates.
       * @param constraints Resolved identifications.
       * @param fixed Whether the fixed degrees of freedom are eliminated.
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
       * @param cache Buffers holding the list.
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
       * @param cache Buffers holding the contributions and the right-hand
       * side.
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
       * @param cache Buffers holding the entries and the right-hand side.
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
       * @param cache Buffers holding the entries and the list.
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
       * @param cache Buffers of the entries and of the list.
       * @param fixed Whether degrees of freedom are fixed.
       * @param x Expansions on the masters, or null without identifications.
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
       * @param vec Vector copied.
       * @param res Destination, with the local size of the vector.
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
       * @param cache Pattern last given to PETSc and merged values.
       * @param A Matrix filled.
       * @param expanded Whether the coordinates are expanded on the masters.
       * @param coordRows Row coordinates of the list.
       * @param coordCols Column coordinates of the list.
       * @param order Coordinates sorted by row then column, as indices in the
       * list.
       * @param starts Start in the sorted coordinates of each distinct
       * coordinate, and the list size.
       * @param values Values of the list.
       */
      template <class Indices, class Scalars>
      void setMatrix(CacheType& cache, ::Mat A, bool expanded, const Indices& coordRows,
        const Indices& coordCols, const COOOrder<MemorySpace>& order,
        const ::Kokkos::View<Index*, MemorySpace>& starts, const Scalars& values) const
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
          if (!expanded && cache.released)
          {
            ::Kokkos::deep_copy(m_execution, rows, cache.distinctRows);
            ::Kokkos::deep_copy(m_execution, cols, cache.distinctCols);
          }
          else
          {
            order.visit([&](const auto& o) {
              listDistinct(m_execution, coordRows, coordCols, o, starts, rows, cols);
            });
          }
          // PETSc may modify the coordinates it is given, which are temporaries.
          const auto hostRows =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), rows);
          const auto hostCols =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), cols);
          if (m_distribution)
          {
            const auto rowGlobals = m_distribution->getPETScIndices();
            const auto colGlobals = getColumns().getPETScIndices();
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
       * @param b Vector filled.
       * @param rhs Right-hand side, of at least the row count.
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
       * @return Whether any process fixes one.
       */
      bool isFixed(bool local) const
      {
        int in = local ? 1 : 0, out = 0;
        MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_LOR, m_distribution->comm);
        return out != 0;
      }

      /**
       * @brief Gets the ghosted vector of two entries per row summing the
       * fixed degrees of freedom over the processes; collective.
       *
       * The vector is kept in the cache, and created again when the layout
       * of the distribution changes on a process, together with the local
       * degrees of freedom held by several processes: the ghosts and the
       * owned degrees of freedom ghosted by another process.
       * @param cache Cache keeping the vector, the shared degrees of freedom
       * and their blocks in the local form.
       * @return Ghosted vector.
       */
      ::Vec getExchange(CacheType& cache) const
      {
        const auto& dist = *m_distribution;
        const int stale = !cache.exchange || !cache.exchangeLayout->hasLayout(dist);
        int any = 0;
        MPI_Allreduce(&stale, &any, 1, MPI_INT, MPI_LOR, dist.comm);
        if (any)
        {
          const PetscInt owned = dist.getOwnedSize();
          const auto rows = dist.getPETScIndices();
          std::vector<PetscInt> ghosts;
          std::vector<PetscInt> blocks(m_rows);
          for (size_t i = 0; i < m_rows; ++i)
          {
            const PetscInt g = rows[i];
            if (g >= dist.begin && g < dist.end)
            {
              blocks[i] = g - dist.begin;
            }
            else
            {
              blocks[i] = owned + static_cast<PetscInt>(ghosts.size());
              ghosts.push_back(g);
            }
          }
          PetscErrorCode ierr;
          ::Vec v;
          ierr = VecCreateGhostBlock(dist.comm, 2, 2 * owned, 2 * dist.size,
            static_cast<PetscInt>(ghosts.size()), ghosts.data(), &v);
          assert(ierr == PETSC_SUCCESS);
          cache.exchange.reset(v, [](::Vec x) { VecDestroy(&x); });
          cache.exchangeLayout = std::make_shared<const COODistribution>(dist);

          // The ghosts count on their owners, which finds the owned degrees
          // of freedom ghosted by another process.
          ::Vec form;
          ierr = VecGhostGetLocalForm(v, &form);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecSet(form, PetscScalar(0));
          assert(ierr == PETSC_SUCCESS);
          PetscScalar* array = nullptr;
          ierr = VecGetArray(form, &array);
          assert(ierr == PETSC_SUCCESS);
          for (size_t k = 0; k < ghosts.size(); ++k)
            array[2 * (owned + static_cast<PetscInt>(k))] = PetscScalar(1);
          ierr = VecRestoreArray(form, &array);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGhostUpdateBegin(v, ADD_VALUES, SCATTER_REVERSE);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGhostUpdateEnd(v, ADD_VALUES, SCATTER_REVERSE);
          assert(ierr == PETSC_SUCCESS);
          const PetscScalar* counts = nullptr;
          ierr = VecGetArrayRead(form, &counts);
          assert(ierr == PETSC_SUCCESS);
          cache.exchangeShared.clear();
          cache.exchangeBlocks.clear();
          for (size_t i = 0; i < m_rows; ++i)
          {
            if (blocks[i] >= owned || PetscRealPart(counts[2 * blocks[i]]) > Real(0.5))
            {
              cache.exchangeShared.push_back(static_cast<Index>(i));
              cache.exchangeBlocks.push_back(blocks[i]);
            }
          }
          ierr = VecRestoreArrayRead(form, &counts);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGhostRestoreLocalForm(v, &form);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
        return cache.exchange.get();
      }

      /**
       * @brief Makes the fixed degrees of freedom the same on every process
       * holding them; collective.
       *
       * A degree of freedom is fixed when one process fixes it, to the mean
       * of the values the processes fixing it impose. The fixed count and
       * the sum of the values of each degree of freedom are summed in one
       * exchange, which only changes the degrees of freedom held by several
       * processes.
       * @param cache Cache keeping the exchange.
       * @param[in,out] mask Whether each local degree of freedom is fixed.
       * @param[in,out] values Fixed value of each local degree of freedom.
       */
      void synchronize(CacheType& cache, std::vector<unsigned char>& mask,
        std::vector<Real>& values) const
      {
        ::Vec v = getExchange(cache);
        const auto& shared = cache.exchangeShared;
        const auto& blocks = cache.exchangeBlocks;
        PetscErrorCode ierr;
        ::Vec form;
        ierr = VecGhostGetLocalForm(v, &form);
        assert(ierr == PETSC_SUCCESS);
        PetscScalar* array = nullptr;
        ierr = VecGetArray(form, &array);
        assert(ierr == PETSC_SUCCESS);
        for (size_t j = 0; j < shared.size(); ++j)
        {
          const Index i = shared[j];
          const PetscInt k = 2 * blocks[j];
          array[k] = mask[i] ? PetscScalar(1) : PetscScalar(0);
          array[k + 1] = mask[i] ? PetscScalar(values[i]) : PetscScalar(0);
        }
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
        const PetscScalar* sums = nullptr;
        ierr = VecGetArrayRead(form, &sums);
        assert(ierr == PETSC_SUCCESS);
        for (size_t j = 0; j < shared.size(); ++j)
        {
          const Index i = shared[j];
          const PetscInt k = 2 * blocks[j];
          const Real c = PetscRealPart(sums[k]);
          mask[i] = c > Real(0.5);
          values[i] = mask[i] ? PetscRealPart(sums[k + 1]) / c : Real(0);
        }
        ierr = VecRestoreArrayRead(form, &sums);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostRestoreLocalForm(v, &form);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
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
       * @param cache Cache keeping the exchange of the fixed degrees of
       * freedom.
       * @return Resolved identifications, in the local numbering.
       */
      typename ConstraintsType::Resolved resolveDistributed(
        [[maybe_unused]] CacheType& cache)
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
          synchronize(cache, mask, values);
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
       * @param cache Buffers holding the right-hand side, and the vector it
       * was last preallocated for.
       * @param b Distributed vector filled.
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
          std::vector<PetscInt> globals = m_distribution->getPETScIndices();
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

}

#endif
