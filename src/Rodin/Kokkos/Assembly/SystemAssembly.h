/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_SYSTEMASSEMBLY_H
#define RODIN_KOKKOS_ASSEMBLY_SYSTEMASSEMBLY_H

/**
 * @file
 * @brief COO reduction, symbolic reuse and constrained Kokkos system assembly.
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosSparse_CooMatrix.hpp>

#if defined(KOKKOS_ENABLE_SYCL)
#include <sycl/ext/intel/math.hpp>
#endif

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Kokkos/ForwardDecls.h"
#include "Rodin/Kokkos/Math/LinearSystem.h"
#include "Rodin/Kokkos/Utility.h"
#include "Integrator.h"
#include "Constraints.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Linear system assembled in device memory with value constraints
   * eliminated, then copied to Eigen.
   *
   * The system is the one of the host assemblies up to rounding, since the
   * integrators sharing a pattern sum their local matrices in place and
   * the multithreaded host assemblies sum in an order that depends on the
   * thread count. The nonzero entries of the local integrators, then of
   * the global integrators, in rank order, followed by the entries of the
   * preassembled forms make the list of matrix contributions. Rows of fixed
   * degrees of freedom are dropped and become identity rows holding the
   * fixed value, and the columns of fixed degrees of freedom move to the
   * right-hand side. The contributions to every
   * coefficient are summed in list order, as
   * Eigen::SparseMatrix::setFromTriplets() sums duplicates, and products
   * are rounded before being summed. The entries stay where the
   * integrators store them, their list order following from their
   * position, and only the compressed matrix and the right-hand side are
   * copied to the host.
   *
   * Identification constraints, from identification Dirichlet conditions
   * and periodic conditions, are resolved on the host as in the host
   * assemblies. Each contribution is then expanded on the device on the
   * masters of its row and column, in rank order, the identified columns
   * lifting their affine values to the right-hand side, and the identified
   * rows become the rows of their constraints.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class SystemAssembly
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

      /// @brief Device array.
      template <class T>
      using Array = ::Kokkos::View<T*, MemorySpace>;

      /// @brief Matrix left on the device.
      using BackendMatrix = ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>;

      /// @brief Device bilinear integrator.
      using BilinearType = ::Rodin::Kokkos::Variational::BilinearFormIntegrator<Backend>;

      /// @brief Device linear integrator.
      using LinearType = ::Rodin::Kokkos::Variational::LinearFormIntegrator<Backend>;

      /// @brief Entries written by one source, contiguous in storage.
      struct Segment
      {
        /// @brief Position of the first entry.
        size_t offset;

        /// @brief Number of entries.
        size_t count;

        /// @brief Spacing of the entries, as returned by getSpacing().
        size_t spacing;

        /// @brief Number of its first polytope.
        size_t first;

        /**
         * @brief Whether the zero entries are kept, as for a preassembled
         * operator.
         */
        bool dense;
      };

      /**
       * @brief Rank of the entries in the list of contributions, by position.
       *
       * The entries following the last segment are ranked by position. The
       * polytopes are numbered in rank order, the segments of spacing one and
       * the entries following them counting as one polytope. Entries of one
       * polytope are stored in rank order, so that sorting keys with the
       * number of the polytope as lowest bits leaves equal keys in rank order
       * when the sort is stable.
       */
      struct Ranks
      {
        /// @brief Segments, by increasing position.
        Array<Segment> segments;

        /// @brief Whether a segment keeps its zero entries.
        bool dense = false;

        /// @brief Number of polytopes.
        size_t polytopes = 1;

        /// @brief Bits of the sort keys holding the polytope.
        unsigned shift = 0;

        /**
         * @brief Finds the segment holding a position.
         * @param e Position of an entry in storage.
         * @return Index of the segment holding @p e, or the number of
         * segments when @p e follows the last one.
         */
        KOKKOS_INLINE_FUNCTION size_t find(size_t e) const
        {
          size_t k = 0;
          while (k < segments.extent(0) && e >= segments(k).offset + segments(k).count)
            ++k;
          return k;
        }

        /**
         * @brief Rank of the entry stored at a position.
         * @param e Position of the entry in storage.
         * @return Rank of the entry in the list of contributions.
         */
        KOKKOS_INLINE_FUNCTION Index operator()(Index e) const
        {
          const size_t k = find(e);
          if (k == segments.extent(0))
            return e;
          const Segment& s = segments(k);
          return s.offset +
            ::Rodin::Kokkos::Variational::getEntryRank(
              e - s.offset, s.count, s.spacing);
        }

        /**
         * @brief Number of the polytope of the entry stored at a position.
         * @param e Position of the entry in storage.
         * @return Number of its polytope, the last one for the entries
         * following the last segment.
         */
        KOKKOS_INLINE_FUNCTION Index polytope(size_t e) const
        {
          const size_t k = find(e);
          if (k == segments.extent(0))
            return polytopes - 1;
          const Segment& s = segments(k);
          return s.spacing > 1 ? s.first + (e - s.offset) % s.spacing : s.first;
        }

        /**
         * @brief Sort key of the entry stored at a position.
         * @param key Key of the entry, without the polytope.
         * @param e Position of the entry in storage.
         * @return Key shifted with the number of the polytope as lowest
         * bits, or @p key when the polytopes do not fit.
         */
        KOKKOS_INLINE_FUNCTION Index sortKey(Index key, size_t e) const
        {
          return shift > 0 ? (key << shift) | polytope(e) : key;
        }

        /**
         * @brief Places the polytopes in the sort keys when they fit.
         * @param count Number of keys.
         * @returns Ranks with the sort keys set.
         */
        Ranks sorting(size_t count) const
        {
          unsigned bits = 0;
          while ((size_t(1) << bits) < polytopes)
            ++bits;
          Ranks res = *this;
          res.shift = bits > 0 && ((count - 1) >> (64 - bits)) == 0 ? bits : 0;
          return res;
        }

        /**
         * @brief Whether the entry stored at a position is listed.
         * @param values Values of the entries, by position.
         * @param e Position of the entry in storage.
         * @return True when the entry is nonzero or belongs to a segment
         * keeping its zero entries.
         */
        KOKKOS_INLINE_FUNCTION bool kept(const Array<Real>& values, size_t e) const
        {
          if (values(e) != Real(0))
            return true;
          if (!dense)
            return false;
          const size_t k = find(e);
          return k < segments.extent(0) && segments(k).dense;
        }
      };

      /// @brief Matrix contributions.
      struct Triplets
      {
        /// @brief Native coordinate matrix holding the contributions.
        using Matrix = ::KokkosSparse::CooMatrix<Real, Index,
          ::Kokkos::Device<ExecutionSpace, MemorySpace>>;

        /// @brief Coordinate storage; copies borrow its arrays.
        Matrix matrix;

        /// @brief Values retained after coordinate sorting releases the COO graph.
        Array<Real> coefficients;

        /// @brief Columns retained only for prescribed-value lifting.
        Array<Index> eliminatedColumns;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /// @brief Number of entries.
        size_t size = 0;
      };

      /// @brief Copies an existing backend matrix into its contribution segment.
      struct CopyMatrix
      {
        /// @brief Compressed source matrix.
        typename BackendMatrix::View source;

        /// @brief Destination coordinate entries.
        typename Triplets::Matrix target;

        /// @brief Start of the contribution segment.
        size_t offset;

        /// @brief Global row offset.
        Index row;

        /// @brief Global column offset.
        Index col;

        /// @param i Matrix row.
        KOKKOS_INLINE_FUNCTION void operator()(Index i) const
        {
          for (auto k = source.graph.row_map(i); k < source.graph.row_map(i + 1); ++k)
          {
            target.row()(offset + k) = row + i;
            target.col()(offset + k) = col + source.graph.entries(k);
            target.data()(offset + k) = source.values(k);
          }
        }
      };

      /// @brief Listed entries.
      struct Listed
      {
        /// @brief Sort keys.
        Array<Index> keys;

        /// @brief Positions.
        Array<Index> positions;
      };

      /// @brief Numbers of listed matrix contributions.
      struct Counts
      {
        /// @brief Contributions to coefficients.
        size_t coefficients = 0;

        /// @brief Contributions moved to the right-hand side.
        size_t eliminated = 0;

        /**
         * @brief Adds counts.
         * @param other Counts to add.
         * @return Reference to these counts.
         */
        KOKKOS_INLINE_FUNCTION Counts& operator+=(const Counts& other)
        {
          coefficients += other.coefficients;
          eliminated += other.eliminated;
          return *this;
        }
      };

      /// @brief Entries grouped by key.
      struct Groups
      {
        /**
         * @brief Start of every group in the entries, followed by the entry
         * count.
         */
        Array<Index> starts;

        /// @brief Key of every group.
        Array<Index> keys;

        /**
         * @brief Number of groups.
         * @return Number of distinct keys.
         */
        size_t size() const
        {
          return keys.extent(0);
        }
      };

      /**
       * @brief Ordered structural groups retained by an assembly policy.
       *
       * Input keys, positions and contribution ranks are owned snapshots. Their
       * exact equality, including dimensions and key encoding, validates reuse.
       * Numerical zeros are selected after grouping, so they cannot hide a new
       * coupling or change the order of the contributions that remain.
       */
      struct Pattern
      {
        /// @brief Whether a complete plan has been built.
        bool valid = false;

        /// @brief Test-space dimension.
        size_t rows = 0;

        /// @brief Trial-space dimension.
        size_t cols = 0;

        /// @brief Number of rank bits encoded in the keys.
        unsigned shift = 0;

        /// @brief Unsorted structural keys.
        Array<Index> keys;

        /// @brief Unsorted source positions.
        Array<Index> positions;

        /// @brief Contribution rank at each unsorted source position.
        Array<Index> ranks;

        /// @brief Source positions grouped in deterministic contribution order.
        Array<Index> ordered;

        /// @brief Starts and keys of the structural groups.
        Groups groups;
      };

      /**
       * @brief Symbolic plans owned by one assembler, in grouping call order.
       *
       * ## Architecture
       *
       * Each execution regenerates structural coordinates. A plan is reused
       * only after an exact backend comparison of every coordinate, source
       * position and rank. Changed layouts rebuild their groups. Values and
       * constraints are always evaluated again. Copies start empty and never
       * share writable working storage with another assembler.
       */
      struct Symbolic
      {
        /// @brief Independent plans, indexed by grouping call order.
        std::vector<Pattern> patterns;

        Symbolic() = default;

        /// @param other Source whose plans are intentionally not shared.
        Symbolic([[maybe_unused]] const Symbolic& other) {}

        /// @param other Source whose plans are intentionally not shared. @returns This empty owner.
        Symbolic& operator=(const Symbolic& other)
        {
          if (this != &other)
            patterns.clear();
          return *this;
        }

        /// @param other Source whose plans are transferred.
        Symbolic(Symbolic&& other) = default;

        /// @param other Source whose plans are transferred. @returns This owner.
        Symbolic& operator=(Symbolic&& other) = default;
      };

      /// @brief Contribution moved to the right-hand side by a fixed column.
      struct Eliminated
      {
        /// @brief Contributions.
        Array<Real> values;

        /// @brief Column indices.
        Array<Index> cols;

        /// @brief Fixed values.
        Array<Real> fixedValues;

        /**
         * @brief Contribution of the entry stored at a position.
         *
         * The product is rounded before being summed, as on the host: device
         * compilers would otherwise fuse it with the sum.
         *
         * @param e Position of the entry in storage.
         * @return Opposite of the product of the entry and the fixed value
         * of its column.
         */
        KOKKOS_INLINE_FUNCTION Real operator()(Index e) const
        {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
          KOKKOS_IF_ON_DEVICE((return -__dmul_rn(values(e), fixedValues(cols(e)));))
#elif defined(KOKKOS_ENABLE_SYCL)
          KOKKOS_IF_ON_DEVICE((return -::sycl::ext::intel::math::dmul_rn(
            values(e), fixedValues(cols(e)));))
#endif
          return -(values(e) * fixedValues(cols(e)));
        }
      };

      /// @brief Constraints of the system.
      using ConstraintsType = Constraints<Backend>;

      /// @brief Expansions of the degrees of freedom on their masters.
      using Expansions = typename ConstraintsType::Expansions;

      /// @brief Numbers of expanded contributions.
      struct Sizes
      {
        /// @brief Matrix or vector entries.
        size_t entries = 0;

        /// @brief Affine values lifted to the right-hand side.
        size_t lifts = 0;

        /**
         * @brief Adds sizes.
         * @param other Sizes to add.
         * @return Reference to these sizes.
         */
        KOKKOS_INLINE_FUNCTION Sizes& operator+=(const Sizes& other)
        {
          entries += other.entries;
          lifts += other.lifts;
          return *this;
        }
      };

      /// @brief Whether a sorted key starts a group.
      struct GroupHead
      {
        /// @brief Sorted keys.
        Array<Index> keys;

        /**
         * @brief Bits of the keys holding the polytope, ignored when
         * comparing.
         */
        unsigned shift;

        /**
         * @brief Whether a sorted key starts a group.
         * @param i Index of the key.
         * @return True when @p i is the first index or its key differs
         * from the previous one.
         */
        KOKKOS_INLINE_FUNCTION bool operator()(const size_t i) const
        {
          return i == 0 || (keys(i) >> shift) != (keys(i - 1) >> shift);
        }
      };

      /**
       * @brief Lists of a contribution, or of a fixed diagonal entry past
       * the triplets.
       */
      struct Classification
      {
        /// @brief Number of contributions, the diagonal entries following.
        size_t size;

        /// @brief Row indices of the contributions.
        Array<Index> rows;

        /// @brief Column indices of the contributions.
        Array<Index> cols;

        /// @brief Values of the contributions.
        Array<Real> values;

        /// @brief Whether each degree of freedom is fixed.
        Array<unsigned char> fixed;

        /// @brief Ranks of the contributions, telling the kept entries.
        Ranks coefficientRanks;

        /// @brief Whether some degree of freedom is fixed.
        bool constrained;

        /// @brief Whether the coefficients are listed.
        bool listCoefficients;

        /// @brief Whether the eliminated contributions are listed.
        bool listEliminated;

        /// @brief Whether zero-valued potential contributions are listed.
        bool structural = false;

        /**
         * @brief Lists of the entry stored at a position.
         * @param e Position of the contribution, or size plus the row of
         * a fixed diagonal entry.
         * @return One coefficient or one eliminated contribution, or
         * none when the entry is dropped.
         */
        KOKKOS_INLINE_FUNCTION Counts operator()(const size_t e) const
        {
          Counts c;
          if (e >= size)
          {
            c.coefficients = fixed(e - size) != 0;
            return c;
          }
          if (!structural && !coefficientRanks.kept(values, e))
            return c;
          const bool row = constrained && fixed(rows(e));
          const bool col = constrained && fixed(cols(e));
          c.coefficients = listCoefficients && !row && !col;
          c.eliminated = listEliminated && !row && col;
          return c;
        }
      };

      /// @brief Expanded sizes of the matrix contribution of rank k.
      struct MatrixExpansionSizes
      {
        /// @brief Position of the contribution of every rank.
        Array<Index> positions;

        /// @brief Row indices of the contributions.
        Array<Index> rows;

        /// @brief Column indices of the contributions.
        Array<Index> cols;

        /// @brief Values of the contributions.
        Array<Real> values;

        /// @brief Expansions of the degrees of freedom.
        Expansions x;

        /**
         * @brief Expanded sizes of a contribution.
         * @param k Rank of the contribution.
         * @return Number of expanded entries, and of affine values lifted
         * by its column, zero for a zero contribution.
         */
        KOKKOS_INLINE_FUNCTION Sizes operator()(const size_t k) const
        {
          Sizes c;
          const Index e = positions(k);
          if (values(e) == Real(0))
            return c;
          const size_t nr = x.count(rows(e));
          c.entries = nr * x.count(cols(e));
          c.lifts = x.values(cols(e)) != Real(0) ? nr : 0;
          return c;
        }
      };

      /// @brief Expanded size of the vector contribution of rank k.
      struct VectorExpansionSize
      {
        /// @brief Position of the contribution of every rank.
        Array<Index> positions;

        /// @brief Row indices of the contributions.
        Array<Index> indices;

        /// @brief Values of the contributions.
        Array<Real> values;

        /// @brief Expansions of the degrees of freedom.
        Expansions x;

        /**
         * @brief Expanded size of a contribution.
         * @param k Rank of the contribution.
         * @return Number of masters of its row, zero for a zero
         * contribution.
         */
        KOKKOS_INLINE_FUNCTION size_t operator()(const size_t k) const
        {
          const Index e = positions(k);
          return values(e) != Real(0) ? x.count(indices(e)) : 0;
        }
      };

      /// @brief Expanded matrix contributions.
      struct Lifted
      {
        /**
         * @brief Contributions, followed by the rows of the identified
         * degrees of freedom.
         */
        Triplets triplets;

        /// @brief Rows of the lifted affine values.
        Array<Index> rows;

        /// @brief Lifted affine values.
        Array<Real> values;

        /// @brief Number of lifts before every rank, followed by their count.
        Array<Index> starts;
      };

      /**
       * @brief Creates an empty system.
       * @param execution Execution space instance of the mesh.
       * @param rows Number of rows, the size of the test spaces.
       * @param cols Number of columns, the size of the trial spaces.
       * @param symbolic Retained plans, or null for a fresh assembly.
       */
      SystemAssembly(const ExecutionSpace& execution, size_t rows, size_t cols,
        Symbolic* symbolic = nullptr)
        : m_execution(execution),
          m_rows(rows),
          m_cols(cols),
          m_constraints(rows, std::max(rows, cols)),
          m_symbolic(symbolic),
          m_group(0)
      {}

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
       * @brief Adds the entries of a local integrator to the matrix.
       * @param lbfi Local integrator, possibly negated or scaled by a real.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void add(
        Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Index row, Index col)
      {
        const auto [integrator, sign] = getIntegrator<Backend>(lbfi);
        m_matrices.push_back({integrator, nullptr, row, col, sign});
      }

      /**
       * @brief Adds the entries of a global integrator to the matrix.
       * @param gbfi Global integrator, possibly negated.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void add(
        Variational::GlobalBilinearFormIntegratorBase<Real>& gbfi, Index row, Index col)
      {
        const auto [integrator, sign] = getIntegrator<Backend>(gbfi);
        m_matrices.push_back({integrator, nullptr, row, col, sign});
      }

      /**
       * @brief Adds the entries of a preassembled operator to the matrix.
       * @param op Operator, alive until finish() returns.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void add(const ::Rodin::Math::SparseMatrix<Real>& op, Index row, Index col)
      {
        m_matrices.push_back({nullptr, &op, row, col, Real(1)});
      }

      /// @param op Backend operator. @param row Test offset. @param col Trial offset.
      void add(const BackendMatrix& op, Index row, Index col)
      {
        m_matrices.push_back({nullptr, nullptr, row, col, 1, &op});
      }

      /// @param vec Backend coefficients. @param row Test offset.
      void add(const ::Rodin::Kokkos::Math::Vector<Real, Backend>& vec, Index row)
      {
        m_vectors.push_back({nullptr, nullptr, row, 1, m_matrices.size(), &vec});
      }

      /**
       * @brief Adds the entries of a linear integrator to the right-hand side.
       *
       * The affine values lifted by the matrices added before are summed
       * first, as the host assemblies lift them while adding the matrices.
       * @param lfi Integrator, possibly negated or scaled by a real.
       * @param row Offset of its test space.
       * @param sign Factor of the contributions.
       */
      void add(Variational::LinearFormIntegratorBase<Real>& lfi, Index row, Real sign)
      {
        const auto [integrator, s] = getIntegrator<Backend>(lfi, sign);
        m_vectors.push_back({integrator, nullptr, row, s, m_matrices.size()});
      }

      /**
       * @brief Adds a preassembled vector to the right-hand side.
       * @param vec Vector, alive until finish() returns.
       * @param row Offset of its test space.
       */
      void add(const ::Rodin::Math::Vector<Real>& vec, Index row)
      {
        m_vectors.push_back({nullptr, &vec, row, Real(1), m_matrices.size()});
      }

      /**
       * @brief Integrates, eliminates the constraints and builds the system.
       * @param A Matrix, built when not null.
       * @param b Right-hand side, built when not null.
       */
      void finish(::Rodin::Math::SparseMatrix<Real>* A, ::Rodin::Math::Vector<Real>* b)
      {
        build(A, b);
      }

      /**
       * @brief Integrates, eliminates the constraints and builds the system,
       * leaving the matrix on the device.
       * @param A Matrix, in compressed rows.
       * @param b Right-hand side.
       *
       * The arrays of the matrix are reused when their sizes are unchanged.
       * The coefficients are summed in the same order as for finish() into a
       * host matrix, so both matrices are equal.
       */
      void finish(BackendMatrix& A, ::Rodin::Math::Vector<Real>& b)
      {
        build(&A, &b);
      }

      /**
       * @brief Builds an entirely backend-resident system.
       * @param A Compressed row operator. @param b Backend right-hand side.
       */
      void finish(BackendMatrix& A, ::Rodin::Kokkos::Math::Vector<Real, Backend>& b)
      {
        build(&A, &b);
      }

      /// @param A Backend operator, optionally absent. @param b Backend right-hand side, optionally absent.
      void finish(BackendMatrix* A, ::Rodin::Kokkos::Math::Vector<Real, Backend>* b)
      {
        build(A, b);
      }

      /// @param A Backend operator. @param b Absent host right-hand side for a bilinear form.
      void finish(BackendMatrix* A, ::Rodin::Math::Vector<Real>* b)
      {
        build(A, b);
      }

      /// @param A Absent operator for a linear form. @param b Backend right-hand side.
      void finish(std::nullptr_t A, ::Rodin::Kokkos::Math::Vector<Real, Backend>* b)
      {
        build(static_cast<BackendMatrix*>(A), b);
      }

    private:
      /**
       * @brief Builds the system into a host or a device matrix.
       * @param A Matrix, host or device, built when not null.
       * @param b Right-hand side, built when not null.
       */
      template <class Matrix, class Vector>
      void build(Matrix* A, Vector* b)
      {
        m_group = 0;
        constexpr bool rowMajor = std::is_same_v<Matrix, BackendMatrix>;
        typename ConstraintsType::Resolved constraints;
        if (m_constraints.isIdentified())
          constraints = m_constraints.resolve(m_execution);
        const bool expanded = constraints.isExpanded();
        const bool constrained = m_constraints.isConstrained();
        Array<unsigned char> fixed;
        Array<Real> fixedValues;
        if (constrained || expanded)
        {
          fixed = upload(m_constraints.getFixed(), "Rodin fixed DOFs");
          fixedValues = upload(m_constraints.getValues(), "Rodin fixed values");
        }
        // The identity rows come last, ranked by position.
        const size_t diagonal = A && constrained ? std::min(m_rows, m_cols) : 0;
        const bool eliminate = b && constrained;
        const auto coefficientRanks = [&](const Ranks& ranks) {
          return ranks.sorting(m_rows * m_cols);
        };
        Triplets triplets;
        Lifted lifted;
        std::vector<size_t> offsets;
        Listed coefficients, eliminated;
        if (A || eliminate || (b && expanded))
        {
          triplets = list(offsets);
          if (expanded)
          {
            lifted = expand(triplets, constraints);
            triplets = lifted.triplets;
          }
          listContributions(m_execution, rowMajor ? m_cols : m_rows, rowMajor, diagonal,
            triplets, fixed, coefficientRanks(triplets.ranks),
            triplets.ranks.sorting(m_rows), A != nullptr, eliminate, coefficients,
            eliminated, true);
          // Sort keys now encode the coordinates. Retain only the arrays
          // still needed by compression and prescribed-value lifting.
          triplets.coefficients = triplets.matrix.data();
          if (eliminate)
            triplets.eliminatedColumns = triplets.matrix.col();
          triplets.matrix = {};
          lifted.triplets.matrix = {};
        }
        if (A)
          compress(triplets, coefficientRanks(triplets.ranks), coefficients, *A);
        if (b)
        {
          const Array<Real> rhs(
            ::Kokkos::view_alloc(m_execution, std::string("Rodin right-hand side")),
            m_rows);
          std::vector<Index> bounds;
          if (expanded)
            bounds = download(lifted.starts, offsets);
          size_t done = 0;
          const auto lift = [&](size_t matrices) {
            if (!expanded || matrices <= done)
              return;
            const auto range = std::make_pair(
              static_cast<size_t>(bounds[done]), static_cast<size_t>(bounds[matrices]));
            done = matrices;
            if (range.first < range.second)
            {
              add(Array<Index>(::Kokkos::subview(lifted.rows, range)),
                Array<Real>(::Kokkos::subview(lifted.values, range)), Ranks(), nullptr,
                rhs);
            }
          };
          sum(rhs, expanded ? &constraints.expansions : nullptr, lift);
          lift(m_matrices.size());
          if (constrained)
          {
            const auto groups = groupContributions(eliminated,
              triplets.ranks.sorting(m_rows), triplets.coefficients, triplets.size);
            accumulate(m_execution, eliminated.positions, groups,
              Eliminated{triplets.coefficients, triplets.eliminatedColumns, fixedValues},
              rhs);
          }
          if (expanded)
          {
            impose(m_execution, m_rows, upload(constraints.imposed, "Rodin imposed DOFs"),
              fixedValues, rhs);
          }
          else if (constrained)
            impose(m_execution, m_rows, fixed, fixedValues, rhs);
          b->resize(m_rows);
          if constexpr (requires { b->getView(); })
            ::Kokkos::deep_copy(m_execution, b->getView(), rhs);
          else
            ::Kokkos::deep_copy(m_execution, host(b->data(), m_rows), rhs);
        }
        m_execution.fence("Rodin linear system");
        if (m_symbolic)
          m_symbolic->patterns.resize(m_group);
      }

    public:
      /// @brief Adds the offsets of their spaces to the indices of an entry.
      struct IndexShift
      {
        /// @brief Row indices, shifted in place.
        Array<Index> rows;

        /// @brief Column indices, shifted in place for a matrix.
        Array<Index> cols;

        /// @brief Offset of the test space.
        Index row;

        /// @brief Offset of the trial space.
        Index col;

        /// @brief Whether the entries have columns.
        bool matrix;

        /**
         * @brief Shifts the indices of an entry.
         * @param e Position of the entry in storage.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          rows(e) += row;
          if (matrix)
            cols(e) += col;
        }
      };

      /**
       * @brief Adds the offsets of their spaces to the indices of entries.
       * @param execution Execution space of the kernel.
       * @param offset Position of the first entry.
       * @param count Number of entries.
       * @param row Offset of the test space.
       * @param col Offset of the trial space.
       * @param rows Row indices.
       * @param cols Column indices, empty for a vector.
       */
      static void shift(const ExecutionSpace& execution, size_t offset, size_t count,
        Index row, Index col, const Array<Index>& rows, const Array<Index>& cols)
      {
        const bool matrix = cols.extent(0) > 0;
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::shift",
          Range(execution, offset, offset + count),
          IndexShift{rows, cols, row, col, matrix});
      }

      /**
       * @brief Counts the groups of sorted keys and the entries following
       * an entry of higher rank in their group.
       */
      struct GroupCount
      {
        /// @brief Tells the keys starting a group.
        GroupHead head;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /// @brief Positions of the entries, in key order.
        Array<Index> positions;

        /**
         * @brief Counts a sorted key.
         * @param i Index of the key.
         * @param d Number of groups, incremented when @p i starts one.
         * @param u Number of entries out of rank order, incremented when
         * the entry of @p i ranks before the previous one of its group.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t i, size_t& d, size_t& u) const
        {
          if (head(i))
            ++d;
          else
            u += ranks(positions(i)) < ranks(positions(i - 1));
        }
      };

      /// @brief Writes the start and the key of every group of sorted keys.
      struct GroupStarts
      {
        /// @brief Tells the keys starting a group.
        GroupHead head;

        /// @brief Receives the start of every group.
        Array<Index> starts;

        /// @brief Receives the key of every group, without the polytope.
        Array<Index> groupKeys;

        /// @brief Sorted keys.
        Array<Index> keys;

        /// @brief Bits of the keys holding the polytope.
        unsigned shift;

        /**
         * @brief Scans a sorted key.
         * @param i Index of the key.
         * @param update Number of groups started before @p i.
         * @param final Whether the pass writes the groups.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t i, size_t& update, const bool final) const
        {
          if (head(i))
          {
            if (final)
            {
              starts(update) = i;
              groupKeys(update) = keys(i) >> shift;
            }
            ++update;
          }
        }
      };

      /// @brief Sorts the positions of a group in rank order, by insertion.
      struct GroupOrder
      {
        /// @brief Start of every group, followed by the entry count.
        Array<Index> starts;

        /// @brief Positions of the entries, sorted in place.
        Array<Index> positions;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /**
         * @brief Sorts a group.
         * @param u Index of the group.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          for (Index i = starts(u) + 1; i < starts(u + 1); ++i)
          {
            const Index p = positions(i);
            const Index r = ranks(p);
            Index j = i;
            for (; j > starts(u) && ranks(positions(j - 1)) > r; --j)
              positions(j) = positions(j - 1);
            positions(j) = p;
          }
        }
      };

      /**
       * @brief Counts the entries of a retained symbolic pattern that differ
       * from the current listing.
       */
      struct SymbolicChanges
      {
        /// @brief Current sort keys.
        Array<Index> keys;

        /// @brief Current positions.
        Array<Index> positions;

        /// @brief Current ranks of the entries.
        Ranks ranks;

        /// @brief Retained sort keys.
        Array<Index> previousKeys;

        /// @brief Retained positions.
        Array<Index> previousPositions;

        /// @brief Retained ranks, in key order.
        Array<Index> previousRanks;

        /**
         * @brief Compares an entry.
         * @param i Index of the entry.
         * @param changed Number of differing entries, incremented.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i, size_t& changed) const
        {
          changed += keys(i) != previousKeys(i) ||
            positions(i) != previousPositions(i) ||
            ranks(positions(i)) != previousRanks(i);
        }
      };

      /// @brief Retains the ranks of the entries in key order.
      struct SymbolicRanks
      {
        /// @brief Receives the ranks.
        Array<Index> saved;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /// @brief Positions of the entries, in key order.
        Array<Index> positions;

        /**
         * @brief Retains the rank of an entry.
         * @param i Index of the entry in key order.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          saved(i) = ranks(positions(i));
        }
      };

      /// @brief Tells the entries kept in the graph presented to callers.
      struct ActiveEntries
      {
        /// @brief Positions of the entries, in group order.
        Array<Index> positions;

        /// @brief Values of the entries, by position.
        Array<Real> values;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /// @brief Positions from this bound on denote identity entries.
        size_t size;

        /**
         * @brief Tells whether an entry is kept.
         * @param k Index of the entry in group order.
         * @return True for an identity entry or a kept value.
         */
        KOKKOS_INLINE_FUNCTION bool operator()(const Index k) const
        {
          const Index e = positions(k);
          return static_cast<size_t>(e) >= size || ranks.kept(values, e);
        }
      };

      /// @brief Counts the kept entries and the groups holding one.
      struct ActiveGroups
      {
        /// @brief Start of every group, followed by the entry count.
        Array<Index> starts;

        /// @brief Tells the kept entries.
        ActiveEntries kept;

        /**
         * @brief Counts a group.
         * @param u Index of the group.
         * @param sum Kept entries in coefficients, and groups in eliminated.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u, Counts& sum) const
        {
          size_t entries = 0;
          for (Index k = starts(u); k < starts(u + 1); ++k)
            entries += kept(k);
          sum.coefficients += entries;
          sum.eliminated += entries != 0;
        }
      };

      /// @brief Writes the groups holding kept entries, with these entries.
      struct GroupSelection
      {
        /// @brief Start of every group, followed by the entry count.
        Array<Index> starts;

        /// @brief Key of every group.
        Array<Index> keys;

        /// @brief Tells the kept entries.
        ActiveEntries kept;

        /// @brief Receives the positions of the kept entries.
        Array<Index> order;

        /// @brief Receives the start of every selected group.
        Array<Index> selectedStarts;

        /// @brief Receives the key of every selected group.
        Array<Index> selectedKeys;

        /**
         * @brief Scans a group.
         * @param u Index of the group.
         * @param sum Kept entries and selected groups before @p u.
         * @param final Whether the pass writes the selection.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t u, Counts& sum, const bool final) const
        {
          size_t entries = 0;
          for (Index k = starts(u); k < starts(u + 1); ++k)
          {
            if (kept(k))
            {
              if (final)
                order(sum.coefficients + entries) = kept.positions(k);
              ++entries;
            }
          }
          if (entries && final)
          {
            selectedStarts(sum.eliminated) = sum.coefficients;
            selectedKeys(sum.eliminated) = keys(u);
          }
          sum.coefficients += entries;
          sum.eliminated += entries != 0;
        }
      };

      /// @brief Writes the sort key and the position of vector entries.
      struct VectorKeys
      {
        /// @brief Rows of the entries.
        Array<Index> indices;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /// @brief Receives the sort keys.
        Array<Index> keys;

        /// @brief Receives the positions.
        Array<Index> positions;

        /**
         * @brief Lists an entry.
         * @param e Position of the entry.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          keys(e) = ranks.sortKey(indices(e), e);
          positions(e) = e;
        }
      };

      /**
       * @brief Sorts keys with the positions of their entries.
       * @param execution Execution space of the kernels.
       * @param keys Sort keys, sorted on return.
       * @param positions Positions of the entries, permuted as the keys on
       * return, the positions of equal keys being in rank order.
       * @param ranks Ranks of the entries, with the layout of the sort keys.
       * @returns Groups of equal keys.
       *
       * The order among equal keys does not depend on the sort being stable.
       */
      static Groups group(const ExecutionSpace& execution, const Array<Index>& keys,
        const Array<Index>& positions, const Ranks& ranks)
      {
        const size_t count = keys.extent(0);
        if (count > 0)
          ::Kokkos::Experimental::sort_by_key(execution, keys, positions);
        const unsigned shift = ranks.shift;
        const GroupHead head{keys, shift};
        const auto policy = Range(execution, 0, count);
        size_t distinct = 0;
        size_t unordered = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::groups", policy,
          GroupCount{head, ranks, positions}, ::Kokkos::Sum<size_t>(distinct),
          ::Kokkos::Sum<size_t>(unordered));
        Groups res;
        res.starts =
          Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                         std::string("Rodin group starts")),
            distinct + 1);
        res.keys =
          Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                         std::string("Rodin group keys")),
            distinct);
        const auto starts = res.starts;
        const auto groupKeys = res.keys;
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::starts", policy,
          GroupStarts{head, starts, groupKeys, keys, shift});
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(starts, distinct), count);
        if (unordered == 0)
          return res;

        // Groups are small, a few polytopes sharing a pair of degrees of freedom.
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::order",
          Range(execution, 0, distinct), GroupOrder{starts, positions, ranks});
        return res;
      }

      /**
       * @brief Groups structural contributions and selects their current values.
       * @param listed Structural keys and positions; receives the selected order.
       * @param ranks Current source ordering and zero-retention rules.
       * @param values Current numerical contributions.
       * @param size Positions from this bound on denote identity entries.
       * @returns Nonempty numerical groups, in the existing summation order.
       */
      Groups groupContributions(
        Listed& listed, const Ranks& ranks, const Array<Real>& values, size_t size)
      {
        const size_t count = listed.keys.extent(0);
        Groups groups;
        if (m_symbolic)
        {
          if (m_group == m_symbolic->patterns.size())
            m_symbolic->patterns.emplace_back();
          auto& pattern = m_symbolic->patterns[m_group++];
          const auto keys = listed.keys;
          const auto positions = listed.positions;
          bool same = pattern.valid && pattern.rows == m_rows && pattern.cols == m_cols &&
            pattern.shift == ranks.shift && pattern.keys.extent(0) == count;
          if (same)
          {
            size_t changes = 0;
            ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::validateSymbolic",
              Range(m_execution, 0, count),
              SymbolicChanges{
                keys, positions, ranks, pattern.keys, pattern.positions, pattern.ranks},
              changes);
            same = changes == 0;
          }
          if (!same)
          {
            pattern = Pattern();
            pattern.keys = allocate<Index>(count, "Rodin symbolic keys");
            pattern.positions = allocate<Index>(count, "Rodin symbolic positions");
            pattern.ranks = allocate<Index>(count, "Rodin symbolic ranks");
            ::Kokkos::deep_copy(m_execution, pattern.keys, keys);
            ::Kokkos::deep_copy(m_execution, pattern.positions, positions);
            ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::symbolicRanks",
              Range(m_execution, 0, count),
              SymbolicRanks{pattern.ranks, ranks, positions});
            pattern.groups = group(m_execution, keys, positions, ranks);
            pattern.ordered = positions;
            pattern.rows = m_rows;
            pattern.cols = m_cols;
            pattern.shift = ranks.shift;
            pattern.valid = true;
          }
          listed.positions = pattern.ordered;
          groups = pattern.groups;
        }
        else
          groups = group(m_execution, listed.keys, listed.positions, ranks);

        // The graph presented to callers retains the existing zero-dropping
        // semantics. Only the retained symbolic plan contains potential zeros.
        const ActiveEntries kept{listed.positions, values, ranks, size};
        Counts active;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::activeGroups",
          Range(m_execution, 0, groups.size()), ActiveGroups{groups.starts, kept},
          active);
        if (active.coefficients == count)
          return groups;
        Groups selected;
        selected.starts = allocate<Index>(active.eliminated + 1, "Rodin active starts");
        selected.keys = allocate<Index>(active.eliminated, "Rodin active keys");
        const auto order = allocate<Index>(active.coefficients, "Rodin active positions");
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::selectGroups",
          Range(m_execution, 0, groups.size()),
          GroupSelection{
            groups.starts, groups.keys, kept, order, selected.starts, selected.keys});
        ::Kokkos::deep_copy(m_execution,
          ::Kokkos::subview(selected.starts, active.eliminated), active.coefficients);
        listed.positions = order;
        return selected;
      }

      /// @brief Sums the lists of the contributions.
      struct ContributionCount
      {
        /// @brief Lists of each contribution.
        Classification classify;

        /**
         * @brief Adds the lists of a contribution.
         * @param e Position of the contribution.
         * @param sum Numbers of listed contributions.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, Counts& sum) const
        {
          sum += classify(e);
        }
      };

      /**
       * @brief Writes the sort key and the position of every listed
       * contribution.
       */
      struct ContributionList
      {
        /// @brief Lists of each contribution.
        Classification classify;

        /// @brief Number of contributions, the diagonal entries following.
        size_t size;

        /// @brief Row indices of the contributions.
        Array<Index> rows;

        /// @brief Column indices of the contributions.
        Array<Index> cols;

        /// @brief Whether the coefficient keys are row major.
        bool rowMajor;

        /**
         * @brief Number of columns of a row major matrix, of rows of a
         * column major one.
         */
        size_t innerCount;

        /// @brief Receives the sort keys of the coefficients.
        Array<Index> ck;

        /// @brief Receives the positions of the coefficients.
        Array<Index> cp;

        /// @brief Ranks with the layout of the coefficient sort keys.
        Ranks coefficientRanks;

        /// @brief Receives the sort keys of the eliminated contributions.
        Array<Index> ek;

        /// @brief Receives the positions of the eliminated contributions.
        Array<Index> ep;

        /// @brief Ranks with the layout of the row sort keys.
        Ranks eliminatedRanks;

        /**
         * @brief Writes the lists of a contribution.
         * @param e Position of the contribution, or size plus the row of
         * a fixed diagonal entry.
         * @param update Numbers of contributions listed before @p e.
         * @param final Whether the pass writes the lists.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, Counts& update, const bool final) const
        {
          const Counts c = classify(e);
          if (final && c.coefficients)
          {
            const Index row = e < size ? rows(e) : e - size;
            const Index col = e < size ? cols(e) : e - size;
            const Index key =
              rowMajor ? row * innerCount + col : col * innerCount + row;
            ck(update.coefficients) = coefficientRanks.sortKey(key, e);
            cp(update.coefficients) = e;
          }
          if (final && c.eliminated)
          {
            ek(update.eliminated) = eliminatedRanks.sortKey(rows(e), e);
            ep(update.eliminated) = e;
          }
          update += c;
        }
      };

      /**
       * @brief Lists the matrix contributions kept on their coefficients
       * with the identity rows, and the contributions moved to the
       * right-hand side by the columns of the fixed degrees of freedom.
       * @param execution Execution space of the kernels.
       * @param innerCount Number of columns of a row major matrix, of rows
       * of a column major one.
       * @param rowMajor Whether the coefficient keys are row major.
       * @param diagonal Number of identity rows listed.
       * @param triplets Matrix contributions.
       * @param fixed Whether each degree of freedom is fixed, empty when
       * none is.
       * @param coefficientRanks Ranks of the contributions, with the layout
       * of the coefficient sort keys.
       * @param eliminatedRanks Ranks of the contributions, with the layout
       * of the row sort keys.
       * @param listCoefficients Whether the coefficients are listed.
       * @param listEliminated Whether the eliminated contributions are listed.
       * @param coefficients Receives the coefficients, as row or column
       * major sort keys, the diagonal of the fixed row i being at position size + i.
       * @param eliminated Receives the eliminated contributions, as row sort
       * keys.
       * @param structural Whether to list potential zero-valued contributions.
       */
      static void listContributions(const ExecutionSpace& execution, size_t innerCount,
        bool rowMajor, size_t diagonal, const Triplets& triplets,
        const Array<unsigned char>& fixed, const Ranks& coefficientRanks,
        const Ranks& eliminatedRanks, bool listCoefficients, bool listEliminated,
        Listed& coefficients, Listed& eliminated, bool structural = false)
      {
        const size_t size = triplets.size;
        const auto rows = triplets.matrix.row();
        const auto cols = triplets.matrix.col();
        const auto values = triplets.matrix.data();
        const bool constrained = fixed.extent(0) > 0;
        const Classification classify{size, rows, cols, values, fixed, coefficientRanks,
          constrained, listCoefficients, listEliminated, structural};
        const auto policy = Range(execution, 0, size + diagonal);
        Counts count;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::countContributions", policy,
          ContributionCount{classify}, count);
        const auto allocate = [&](Listed& listed, size_t n, const char* label) {
          listed.keys =
            Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                           std::string(label) + " keys"),
              n);
          listed.positions =
            Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                           std::string(label) + " positions"),
              n);
        };
        allocate(coefficients, count.coefficients, "Rodin coefficient");
        allocate(eliminated, count.eliminated, "Rodin eliminated");
        const auto ck = coefficients.keys;
        const auto cp = coefficients.positions;
        const auto ek = eliminated.keys;
        const auto ep = eliminated.positions;
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::listContributions", policy,
          ContributionList{classify, size, rows, cols, rowMajor, innerCount, ck, cp,
            coefficientRanks, ek, ep, eliminatedRanks});
      }

      /**
       * @brief Sums the contributions of a coefficient and writes its inner
       * index.
       */
      struct CoefficientSum
      {
        /// @brief Number of contributions, the diagonal entries following.
        size_t size;

        /// @brief Values of the contributions.
        Array<Real> values;

        /**
         * @brief Positions of the listed contributions, grouped by
         * coefficient in rank order.
         */
        Array<Index> positions;

        /**
         * @brief Start of every coefficient in the positions, followed by
         * their count.
         */
        Array<Index> starts;

        /// @brief Receives the value of every coefficient.
        Array<Real> data;

        /// @brief Receives the inner index of every coefficient.
        Array<::Rodin::Math::SparseIndex> inner;

        /// @brief Sort key of every coefficient, without the polytope.
        Array<Index> keys;

        /**
         * @brief Number of columns of a row major matrix, of rows of a
         * column major one.
         */
        size_t innerCount;

        /**
         * @brief Sums a coefficient.
         * @param u Index of the coefficient in the compressed storage.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          const auto value = [&](const Index e) {
            return e < size ? values(e) : Real(1);
          };
          Real sum = value(positions(starts(u)));
          for (Index i = starts(u) + 1; i < starts(u + 1); ++i)
            sum += value(positions(i));
          data(u) = sum;
          inner(u) = static_cast<::Rodin::Math::SparseIndex>(keys(u) % innerCount);
        }
      };

      /**
       * @brief Starts the outer indices following the one of the previous
       * coefficient, up to its own, at a coefficient.
       */
      struct OuterStarts
      {
        /// @brief Number of coefficients.
        size_t nnz;

        /// @brief Sort key of every coefficient, without the polytope.
        Array<Index> keys;

        /**
         * @brief Number of columns of a row major matrix, of rows of a
         * column major one.
         */
        size_t innerCount;

        /**
         * @brief Number of rows of a row major matrix, of columns of a
         * column major one.
         */
        size_t outerCount;

        /// @brief Receives the start of every outer index.
        Array<::Rodin::Math::SparseIndex> outer;

        /**
         * @brief Writes the outer starts ending at a coefficient.
         * @param u Index of the coefficient, or the number of coefficients
         * for the starts following the last one.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          const Index end = u < nnz ? keys(u) / innerCount : outerCount;
          const Index begin = u > 0 ? keys(u - 1) / innerCount + 1 : 0;
          for (Index j = begin; j <= end; ++j)
            outer(j) = static_cast<::Rodin::Math::SparseIndex>(u);
        }
      };

      /**
       * @brief Sums the contributions of every coefficient into the
       * compressed storage.
       * @param execution Execution space of the kernels.
       * @param innerCount Number of columns of a row major matrix, of rows
       * of a column major one.
       * @param outerCount Number of rows of a row major matrix, of columns
       * of a column major one.
       * @param triplets Matrix contributions.
       * @param positions Positions of the listed contributions, grouped by
       * coefficient in rank order.
       * @param groups Groups of the listed contributions, one by
       * coefficient.
       * @param outer Receives the start of every outer index.
       * @param inner Receives the inner index of every coefficient.
       * @param data Receives the value of every coefficient.
       */
      static void sumCoefficients(const ExecutionSpace& execution, size_t innerCount,
        size_t outerCount, const Triplets& triplets, const Array<Index>& positions,
        const Groups& groups, const Array<::Rodin::Math::SparseIndex>& outer,
        const Array<::Rodin::Math::SparseIndex>& inner, const Array<Real>& data)
      {
        const size_t size = triplets.size;
        const auto values = triplets.coefficients;
        const auto starts = groups.starts;
        const auto keys = groups.keys;
        const size_t nnz = data.extent(0);
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::sumCoefficients",
          Range(execution, 0, nnz),
          CoefficientSum{size, values, positions, starts, data, inner, keys, innerCount});
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::outer",
          Range(execution, 0, nnz + 1),
          OuterStarts{nnz, keys, innerCount, outerCount, outer});
      }

      /// @brief Counts the nonzero entries of a vector.
      struct NonzeroCount
      {
        /// @brief Values of the entries.
        Array<Real> values;

        /**
         * @brief Counts an entry when it is nonzero.
         * @param e Position of the entry.
         * @param sum Number of nonzero entries.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e, size_t& sum) const
        {
          sum += values(e) != Real(0);
        }
      };

      /**
       * @brief Writes the sort key and the position of every nonzero entry
       * of a vector.
       */
      struct EntryList
      {
        /// @brief Values of the entries.
        Array<Real> values;

        /// @brief Receives the sort keys of the nonzero entries.
        Array<Index> k;

        /// @brief Receives the positions of the nonzero entries.
        Array<Index> p;

        /// @brief Ranks with the layout of the sort keys.
        Ranks ranks;

        /// @brief Indices of the entries.
        Array<Index> indices;

        /**
         * @brief Lists an entry when it is nonzero.
         * @param e Position of the entry.
         * @param update Number of nonzero entries before @p e.
         * @param final Whether the pass writes the list.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t e, size_t& update, const bool final) const
        {
          if (values(e) != Real(0))
          {
            if (final)
            {
              k(update) = ranks.sortKey(indices(e), e);
              p(update) = e;
            }
            ++update;
          }
        }
      };

      /**
       * @brief Lists the nonzero entries of a vector.
       * @param execution Execution space of the kernels.
       * @param indices Indices of the entries.
       * @param values Values of the entries.
       * @param ranks Ranks of the entries, with the layout of the sort keys.
       * @param keys Receives the sort keys of their indices.
       * @param positions Receives their positions.
       */
      static void listEntries(const ExecutionSpace& execution,
        const Array<Index>& indices, const Array<Real>& values, const Ranks& ranks,
        Array<Index>& keys, Array<Index>& positions)
      {
        const auto policy = Range(execution, 0, indices.extent(0));
        size_t count = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::countEntries", policy, NonzeroCount{values}, count);
        keys = Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                              std::string("Rodin vector keys")),
          count);
        positions =
          Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                         std::string("Rodin vector positions")),
            count);
        const auto k = keys;
        const auto p = positions;
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::listEntries", policy,
          EntryList{values, k, p, ranks, indices});
      }

      /**
       * @brief Adds the contributions of a group to its entry of the
       * right-hand side, in rank order.
       */
      template <class Contribution>
      struct Accumulation
      {
        /// @brief Row of every group.
        Array<Index> keys;

        /// @brief Right-hand side, updated in place.
        Array<Real> rhs;

        /**
         * @brief Start of every group in the positions, followed by their
         * count.
         */
        Array<Index> starts;

        /// @brief Contribution of an entry, by position.
        Contribution contributions;

        /// @brief Positions of the entries, grouped in rank order.
        Array<Index> positions;

        /**
         * @brief Adds the contributions of a group.
         * @param u Index of the group.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t u) const
        {
          const Index i = keys(u);
          Real x = rhs(i);
          for (Index k = starts(u); k < starts(u + 1); ++k)
            x += contributions(positions(k));
          rhs(i) = x;
        }
      };

      /**
       * @brief Adds grouped contributions to the right-hand side, in rank
       * order.
       * @tparam Contribution Callable giving the contribution of an entry
       * from its position.
       * @param execution Execution space of the kernel.
       * @param positions Positions of the entries, grouped in rank order.
       * @param groups Groups of the entries, one by row.
       * @param contributions Contribution of an entry, by position.
       * @param rhs Right-hand side, updated in place.
       */
      template <class Contribution>
      static void accumulate(const ExecutionSpace& execution,
        const Array<Index>& positions, const Groups& groups,
        const Contribution& contributions, const Array<Real>& rhs)
      {
        const auto starts = groups.starts;
        const auto keys = groups.keys;
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::accumulate",
          Range(execution, 0, groups.size()),
          Accumulation<Contribution>{keys, rhs, starts, contributions, positions});
      }

      /// @brief Adds an entry of a preassembled vector to the right-hand side.
      struct DenseSum
      {
        /// @brief Right-hand side, updated in place.
        Array<Real> rhs;

        /// @brief Offset of the test space.
        Index row;

        /// @brief Entries of the vector.
        Array<Real> values;

        /**
         * @brief Adds an entry.
         * @param i Index of the entry in the vector.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          rhs(row + i) += values(i);
        }
      };

      /**
       * @brief Adds a preassembled vector to the right-hand side.
       * @param execution Execution space of the kernel.
       * @param values Entries of the vector.
       * @param row Offset of its test space.
       * @param rhs Right-hand side, updated in place.
       */
      static void add(const ExecutionSpace& execution, const Array<Real>& values,
        Index row, const Array<Real>& rhs)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::addDense",
          Range(execution, 0, values.extent(0)), DenseSum{rhs, row, values});
      }

      /// @brief Writes the fixed value on a fixed row.
      struct FixedValue
      {
        /// @brief Whether each row is fixed.
        Array<unsigned char> fixed;

        /// @brief Right-hand side, updated in place.
        Array<Real> rhs;

        /// @brief Fixed value of each row.
        Array<Real> values;

        /**
         * @brief Writes the fixed value of a row when it is fixed.
         * @param i Index of the row.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          if (fixed(i))
            rhs(i) = values(i);
        }
      };

      /**
       * @brief Writes the fixed values on the fixed rows.
       * @param execution Execution space of the kernel.
       * @param count Number of rows.
       * @param fixed Whether each row is fixed.
       * @param values Fixed value of each row.
       * @param rhs Right-hand side, updated in place.
       */
      static void impose(const ExecutionSpace& execution, size_t count,
        const Array<unsigned char>& fixed, const Array<Real>& values,
        const Array<Real>& rhs)
      {
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::imposeFixed",
          Range(execution, 0, count), FixedValue{fixed, rhs, values});
      }

      /// @brief Writes the position of an entry at its rank.
      struct RankOrder
      {
        /// @brief Receives the position of the entry of every rank.
        Array<Index> res;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /**
         * @brief Writes the position of an entry.
         * @param e Position of the entry in storage.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t e) const
        {
          res(ranks(e)) = e;
        }
      };

      /**
       * @brief Inverts the ranks of entries.
       * @param execution Execution space of the kernel.
       * @param count Number of entries.
       * @param ranks Ranks of the entries.
       * @returns Position of the entry of every rank.
       */
      static Array<Index> order(
        const ExecutionSpace& execution, size_t count, const Ranks& ranks)
      {
        const Array<Index> res(
          ::Kokkos::view_alloc(
            execution, ::Kokkos::WithoutInitializing, std::string("Rodin rank order")),
          count);
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::rankOrder",
          Range(execution, 0, count), RankOrder{res, ranks});
        return res;
      }

      /// @brief Sums the expanded sizes of the matrix contributions.
      struct MatrixExpansionCount
      {
        /// @brief Expanded sizes of each contribution.
        MatrixExpansionSizes sizes;

        /**
         * @brief Adds the expanded sizes of a contribution.
         * @param k Rank of the contribution.
         * @param sum Total expanded sizes.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t k, Sizes& sum) const
        {
          sum += sizes(k);
        }
      };

      /**
       * @brief Expands the matrix contribution of rank k and lifts its
       * affine values.
       */
      struct MatrixExpansion
      {
        /// @brief Expanded sizes of each contribution.
        MatrixExpansionSizes sizes;

        /// @brief Receives the number of lifts before every rank.
        Array<Index> ls;

        /// @brief Position of the contribution of every rank.
        Array<Index> positions;

        /// @brief Values of the contributions.
        Array<Real> values;

        /// @brief Expansions of the degrees of freedom.
        Expansions x;

        /// @brief Row indices of the contributions.
        Array<Index> rows;

        /// @brief Column indices of the contributions.
        Array<Index> cols;

        /// @brief Receives the rows of the expanded contributions.
        Array<Index> er;

        /// @brief Receives the columns of the expanded contributions.
        Array<Index> ec;

        /// @brief Receives the values of the expanded contributions.
        Array<Real> ev;

        /// @brief Receives the rows of the lifted affine values.
        Array<Index> lr;

        /// @brief Receives the lifted affine values.
        Array<Real> lv;

        /**
         * @brief Expands a contribution.
         * @param k Rank of the contribution.
         * @param update Expanded sizes of the contributions of lower rank.
         * @param final Whether the pass writes the expansion.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t k, Sizes& update, const bool final) const
        {
          const Sizes c = sizes(k);
          if (final)
          {
            ls(k) = update.lifts;
            if (c.entries > 0)
            {
              const Index e = positions(k);
              const Real v = values(e);
              const Index r0 = x.starts(rows(e));
              const Index c0 = x.starts(cols(e));
              const size_t nr = x.count(rows(e));
              const size_t nc = x.count(cols(e));
              const Real d = x.values(cols(e));
              for (size_t ri = 0; ri < nr; ++ri)
              {
                const Real rv = x.coefficients(r0 + ri) * v;
                for (size_t ci = 0; ci < nc; ++ci)
                {
                  const size_t q = update.entries + ri * nc + ci;
                  er(q) = x.indices(r0 + ri);
                  ec(q) = x.indices(c0 + ci);
                  ev(q) = rv * x.coefficients(c0 + ci);
                }
                if (c.lifts > 0)
                {
                  lr(update.lifts + ri) = x.indices(r0 + ri);
                  lv(update.lifts + ri) = -(rv * d);
                }
              }
            }
          }
          update += c;
        }
      };

      /**
       * @brief Expands the matrix contributions on the masters of their rows
       * and columns, in rank order.
       * @param execution Execution space of the kernels.
       * @param triplets Matrix contributions.
       * @param x Expansions of the degrees of freedom.
       * @param extra Number of entries reserved after the expanded ones.
       * @param res Receives the contributions, without their ranks, and the
       * lifted affine values.
       *
       * The nonzero contribution @f$ a @f$ of row @f$ i @f$ and column
       * @f$ j @f$ becomes the contributions @f$ (c_r a) c_s @f$ of the
       * masters @f$ r @f$ of @f$ i @f$ and @f$ s @f$ of @f$ j @f$, and lifts
       * @f$ -(c_r a) d_j @f$ to row @f$ r @f$ when @f$ j @f$ has a nonzero
       * affine value @f$ d_j @f$, as the host assemblies do.
       */
      static void expand(const ExecutionSpace& execution, const Triplets& triplets,
        const Expansions& x, size_t extra, Lifted& res)
      {
        const size_t size = triplets.size;
        const auto rows = triplets.matrix.row();
        const auto cols = triplets.matrix.col();
        const auto values = triplets.matrix.data();
        const auto positions = order(execution, size, triplets.ranks);
        const MatrixExpansionSizes sizes{positions, rows, cols, values, x};
        const auto policy = Range(execution, 0, size);
        Sizes total;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::countExpansions", policy,
          MatrixExpansionCount{sizes}, total);
        const auto allocate = [&](auto& array, size_t n, const char* label) {
          using ArrayType = std::decay_t<decltype(array)>;
          array = ArrayType(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                              std::string(label)),
            n);
        };
        res.triplets.size = total.entries + extra;
        Array<Index> rowsOut, colsOut;
        Array<Real> valuesOut;
        allocate(rowsOut, res.triplets.size, "Rodin expanded rows");
        allocate(colsOut, res.triplets.size, "Rodin expanded columns");
        allocate(valuesOut, res.triplets.size, "Rodin expanded values");
        res.triplets.matrix = typename Triplets::Matrix(triplets.matrix.numRows(),
          triplets.matrix.numCols(), rowsOut, colsOut, valuesOut);
        allocate(res.rows, total.lifts, "Rodin lifted rows");
        allocate(res.values, total.lifts, "Rodin lifted values");
        allocate(res.starts, size + 1, "Rodin lifted starts");
        const auto er = res.triplets.matrix.row();
        const auto ec = res.triplets.matrix.col();
        const auto ev = res.triplets.matrix.data();
        const auto lr = res.rows;
        const auto lv = res.values;
        const auto ls = res.starts;
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::expandMatrix", policy,
          MatrixExpansion{
            sizes, ls, positions, values, x, rows, cols, er, ec, ev, lr, lv});
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(ls, size), total.lifts);
      }

      /// @brief Sums the expanded sizes of the vector contributions.
      struct VectorExpansionCount
      {
        /// @brief Expanded size of each contribution.
        VectorExpansionSize sizes;

        /**
         * @brief Adds the expanded size of a contribution.
         * @param k Rank of the contribution.
         * @param sum Total expanded size.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t k, size_t& sum) const
        {
          sum += sizes(k);
        }
      };

      /// @brief Expands the vector contribution of rank k.
      struct VectorExpansion
      {
        /// @brief Expanded size of each contribution.
        VectorExpansionSize sizes;

        /// @brief Position of the contribution of every rank.
        Array<Index> positions;

        /// @brief Values of the contributions.
        Array<Real> values;

        /// @brief Expansions of the degrees of freedom.
        Expansions x;

        /// @brief Indices of the contributions.
        Array<Index> indices;

        /// @brief Receives the indices of the expanded contributions.
        Array<Index> ei;

        /// @brief Receives the values of the expanded contributions.
        Array<Real> ev;

        /**
         * @brief Expands a contribution.
         * @param k Rank of the contribution.
         * @param update Expanded size of the contributions of lower rank.
         * @param final Whether the pass writes the expansion.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          const size_t k, size_t& update, const bool final) const
        {
          const size_t n = sizes(k);
          if (final && n > 0)
          {
            const Index e = positions(k);
            const Real v = values(e);
            const Index r0 = x.starts(indices(e));
            for (size_t ri = 0; ri < n; ++ri)
            {
              ei(update + ri) = x.indices(r0 + ri);
              ev(update + ri) = x.coefficients(r0 + ri) * v;
            }
          }
          update += n;
        }
      };

      /**
       * @brief Expands the vector contributions on the masters of their
       * indices, in rank order.
       * @param execution Execution space of the kernels.
       * @param indices Indices of the contributions.
       * @param values Values of the contributions.
       * @param ranks Ranks of the contributions.
       * @param x Expansions of the degrees of freedom.
       * @param resIndices Receives the indices of the expanded
       * contributions, by rank.
       * @param resValues Receives the expanded contributions, by rank.
       *
       * The nonzero contribution @f$ a @f$ of index @f$ i @f$ becomes the
       * contributions @f$ c_r a @f$ of the masters @f$ r @f$ of @f$ i @f$.
       */
      static void expand(const ExecutionSpace& execution, const Array<Index>& indices,
        const Array<Real>& values, const Ranks& ranks, const Expansions& x,
        Array<Index>& resIndices, Array<Real>& resValues)
      {
        const size_t size = indices.extent(0);
        const auto positions = order(execution, size, ranks);
        const VectorExpansionSize sizes{positions, indices, values, x};
        const auto policy = Range(execution, 0, size);
        size_t total = 0;
        ::Kokkos::parallel_reduce("Rodin::Assembly::Kokkos::countVectorExpansions",
          policy, VectorExpansionCount{sizes}, total);
        resIndices =
          Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                         std::string("Rodin expanded indices")),
            total);
        resValues =
          Array<Real>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                        std::string("Rodin expanded contributions")),
            total);
        const auto ei = resIndices;
        const auto ev = resValues;
        ::Kokkos::parallel_scan("Rodin::Assembly::Kokkos::expandVector", policy,
          VectorExpansion{sizes, positions, values, x, indices, ei, ev});
      }

      /// @brief Writes the index of an entry of a preassembled vector.
      struct Sequence
      {
        /// @brief Receives the indices.
        Array<Index> res;

        /// @brief Offset of the space of the vector.
        Index first;

        /**
         * @brief Writes an index.
         * @param i Index of the entry in the vector.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
        {
          res(i) = first + i;
        }
      };

      /**
       * @brief Indices of a preassembled vector in the system.
       * @param execution Execution space of the kernel.
       * @param count Size of the vector. @param first Offset of its space.
       * @return Consecutive indices starting at @p first.
       */
      static Array<Index> sequence(
        const ExecutionSpace& execution, size_t count, Index first)
      {
        const Array<Index> res(
          ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
            std::string("Rodin preassembled indices")),
          count);
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::sequence",
          Range(execution, 0, count), Sequence{res, first});
        return res;
      }

    private:
      struct Matrix
      {
        BilinearType* integrator;
        const ::Rodin::Math::SparseMatrix<Real>* op;
        Index row;
        Index col;
        Real sign;
        const BackendMatrix* backend = nullptr;
      };

      struct Vector
      {
        LinearType* integrator;
        const ::Rodin::Math::Vector<Real>* vec;
        Index row;
        Real sign;
        // Number of matrices added before.
        size_t matrices;
        const ::Rodin::Kokkos::Math::Vector<Real, Backend>* backend = nullptr;
      };

      /**
       * @brief Wraps a host buffer in an unmanaged view.
       * @param data Host buffer.
       * @param size Number of entries of @p data.
       * @return Unmanaged host view of the buffer.
       */
      template <class T>
      static ::Kokkos::View<T*, ::Kokkos::HostSpace,
        ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>
      host(T* data, size_t size)
      {
        return ::Kokkos::View<T*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(data, size);
      }

      /**
       * @brief Allocates an uninitialized device array.
       * @param size Number of entries.
       * @param label Label of the array.
       * @return Uninitialized device array.
       */
      template <class T>
      Array<T> allocate(size_t size, const char* label) const
      {
        return Array<T>(::Kokkos::view_alloc(
                          m_execution, ::Kokkos::WithoutInitializing, std::string(label)),
          size);
      }

      /**
       * @brief Copies a host buffer into a new device array, waiting for the
       * copy to complete.
       * @param data Host buffer.
       * @param size Number of entries of @p data.
       * @param label Label of the device array.
       * @return Device array holding a copy of the buffer.
       */
      template <class T>
      Array<T> upload(const T* data, size_t size, const char* label) const
      {
        const auto device = allocate<T>(size, label);
        ::Kokkos::deep_copy(m_execution, device, host(data, size));
        // The host buffer may be released before the copy.
        m_execution.fence("Rodin upload");
        return device;
      }

      /**
       * @brief Copies a host vector into a new device array, waiting for the
       * copy to complete.
       * @param data Host values.
       * @param label Label of the device array.
       * @return Device array holding a copy of @p data.
       */
      template <class T>
      Array<T> upload(const std::vector<T>& data, const char* label) const
      {
        return upload(data.data(), data.size(), label);
      }

      /**
       * @brief Integrates the local integrators and lists the matrix
       * contributions.
       * @param offsets Set to the rank of the first contribution of every
       * matrix, followed by their count.
       * @return Listed contributions of all the matrices, in the order they
       * were added.
       */
      Triplets list(std::vector<size_t>& offsets) const
      {
        std::vector<Segment> segments;
        std::vector<size_t> counts;
        size_t total = 0;
        offsets.assign(1, 0);
        for (const auto& m : m_matrices)
        {
          // Each integrator owns its entries, listed on the diagonals of the
          // component blocks when the others vanish.
          const size_t count = m.integrator
            ? m.integrator->prepare(true)
            : (m.backend ? m.backend->getNonZeroCount() : m.op->nonZeros());
          counts.push_back(count);
          offsets.push_back(total + count);
          if (count > 0)
          {
            segments.push_back({total, count,
              m.integrator ? m.integrator->getSpacing() : size_t(1), 0, !m.integrator});
          }
          total += count;
        }
        Triplets triplets;
        triplets.size = total;
        triplets.matrix = typename Triplets::Matrix(m_rows, m_cols,
          allocate<Index>(total, "Rodin matrix rows"),
          allocate<Index>(total, "Rodin matrix columns"),
          allocate<Real>(total, "Rodin matrix values"));
        size_t offset = 0;
        ::Rodin::Kokkos::Internal::clearLaunches();
        for (size_t k = 0; k < m_matrices.size(); ++k)
        {
          const auto& m = m_matrices[k];
          const size_t count = counts[k];
          if (m.integrator)
          {
            m.integrator->fill(triplets.matrix.data(), offset, m.sign, false,
              triplets.matrix.row(), triplets.matrix.col());
            ::Rodin::Kokkos::Internal::checkLaunches("A device bilinear integrator");
            if (m.row != 0 || m.col != 0)
            {
              shift(m_execution, offset, count, m.row, m.col, triplets.matrix.row(),
                triplets.matrix.col());
            }
            offset += count;
            continue;
          }
          if (m.backend)
          {
            ::Kokkos::parallel_for("Rodin preassembled matrix",
              Range(m_execution, 0, m.backend->getRowCount()),
              CopyMatrix{m.backend->getView(), triplets.matrix, offset, m.row, m.col});
            offset += count;
            continue;
          }
          std::vector<Index> rows, cols;
          std::vector<Real> values;
          rows.reserve(count);
          cols.reserve(count);
          values.reserve(count);
          for (int j = 0; j < m.op->outerSize(); ++j)
          {
            for (typename ::Rodin::Math::SparseMatrix<Real>::InnerIterator it(*m.op, j);
                 it; ++it)
            {
              rows.push_back(m.row + it.row());
              cols.push_back(m.col + it.col());
              values.push_back(it.value());
            }
          }
          const auto range = std::make_pair(offset, offset + count);
          ::Kokkos::deep_copy(m_execution,
            ::Kokkos::subview(triplets.matrix.row(), range),
            host(rows.data(), rows.size()));
          ::Kokkos::deep_copy(m_execution,
            ::Kokkos::subview(triplets.matrix.col(), range),
            host(cols.data(), cols.size()));
          ::Kokkos::deep_copy(m_execution,
            ::Kokkos::subview(triplets.matrix.data(), range),
            host(values.data(), values.size()));
          m_execution.fence("Rodin preassembled operator");
          offset += count;
        }
        triplets.ranks = rank(segments, "Rodin matrix segments");
        return triplets;
      }

      /**
       * @brief Numbers the polytopes of segments and uploads them.
       * @param segments Segments of contributions, numbered consecutively.
       * @param label Label of the device array of segments.
       * @return Ranks over the uploaded segments.
       */
      Ranks rank(std::vector<Segment> segments, const char* label) const
      {
        Ranks ranks;
        ranks.polytopes = 0;
        for (auto& s : segments)
        {
          s.first = ranks.polytopes;
          ranks.polytopes += s.spacing;
          ranks.dense = ranks.dense || s.dense;
        }
        // The entries following the segments.
        ++ranks.polytopes;
        ranks.segments = upload(segments, label);
        return ranks;
      }

      /**
       * @brief Sums the listed coefficients and copies the compressed matrix.
       * @param triplets Matrix contributions.
       * @param ranks Ranks of the contributions, with the layout of the
       * coefficient sort keys.
       * @param listed Sort keys and positions of the contributions to the
       * coefficients.
       * @param A Host matrix, resized and filled in compressed columns.
       */
      void compress(const Triplets& triplets, const Ranks& ranks, Listed& listed,
        ::Rodin::Math::SparseMatrix<Real>& A)
      {
        A.resize(m_rows, m_cols);
        if (m_rows == 0 || m_cols == 0)
          return;
        const auto groups =
          groupContributions(listed, ranks, triplets.coefficients, triplets.size);
        const size_t nnz = groups.size();
        if (nnz >
          static_cast<size_t>(std::numeric_limits<::Rodin::Math::SparseIndex>::max()))
        {
          Alert::Exception() << "The matrix is too large for the sparse index type."
                             << Alert::Raise;
        }
        const auto outer =
          allocate<::Rodin::Math::SparseIndex>(m_cols + 1, "Rodin matrix outer");
        const auto inner =
          allocate<::Rodin::Math::SparseIndex>(nnz, "Rodin matrix inner");
        const auto data = allocate<Real>(nnz, "Rodin matrix data");
        sumCoefficients(m_execution, m_rows, m_cols, triplets, listed.positions, groups,
          outer, inner, data);
        A.resizeNonZeros(static_cast<Eigen::Index>(nnz));
        ::Kokkos::deep_copy(m_execution, host(A.outerIndexPtr(), m_cols + 1), outer);
        ::Kokkos::deep_copy(m_execution, host(A.innerIndexPtr(), nnz), inner);
        ::Kokkos::deep_copy(m_execution, host(A.valuePtr(), nnz), data);
      }

      /**
       * @brief Sums the listed coefficients into compressed rows on the
       * device.
       * @param triplets Matrix contributions.
       * @param ranks Ranks of the contributions, with the layout of the
       * coefficient sort keys.
       * @param listed Sort keys and positions of the contributions to the
       * coefficients.
       * @param A Device matrix, whose arrays are reused when their sizes are
       * unchanged.
       */
      void compress(
        const Triplets& triplets, const Ranks& ranks, Listed& listed, BackendMatrix& A)
      {
        if (m_rows == 0 || m_cols == 0)
        {
          A = BackendMatrix(m_rows, m_cols, {}, {}, {});
          return;
        }
        const auto groups =
          groupContributions(listed, ranks, triplets.coefficients, triplets.size);
        const size_t nnz = groups.size();
        if (nnz >
          static_cast<size_t>(std::numeric_limits<::Rodin::Math::SparseIndex>::max()))
        {
          Alert::Exception() << "The matrix is too large for the sparse index type."
                             << Alert::Raise;
        }
        // A sequence of systems with one pattern keeps its arrays.
        const auto reuse = [this](const auto& a, size_t size, const char* label) {
          using T = typename std::decay_t<decltype(a)>::non_const_value_type;
          return a.extent(0) == size ? Array<T>(a)
                                     : this->template allocate<T>(size, label);
        };
        const auto starts =
          reuse(A.getRowStarts(), m_rows + 1, "Rodin matrix row starts");
        const auto columns = reuse(A.getColumns(), nnz, "Rodin matrix columns");
        const auto values = reuse(A.getValues(), nnz, "Rodin matrix values");
        sumCoefficients(m_execution, m_cols, m_rows, triplets, listed.positions, groups,
          starts, columns, values);
        A = BackendMatrix(m_rows, m_cols, starts, columns, values);
      }

      /**
       * @brief Expands the matrix contributions and appends the rows of the
       * identified degrees of freedom, as the host assemblies do.
       * @param triplets Matrix contributions.
       * @param constraints Resolved constraints, with the expansions on the
       * masters and the rows of the identified degrees of freedom.
       * @return Expanded contributions and the lifted affine values.
       */
      Lifted expand(const Triplets& triplets,
        const typename ConstraintsType::Resolved& constraints) const
      {
        Lifted res;
        const size_t extra = constraints.rows.size();
        expand(m_execution, triplets, constraints.expansions, extra, res);
        const size_t count = res.triplets.size - extra;
        const auto range = std::make_pair(count, res.triplets.size);
        ::Kokkos::deep_copy(m_execution,
          ::Kokkos::subview(res.triplets.matrix.row(), range),
          host(constraints.rows.data(), extra));
        ::Kokkos::deep_copy(m_execution,
          ::Kokkos::subview(res.triplets.matrix.col(), range),
          host(constraints.cols.data(), extra));
        ::Kokkos::deep_copy(m_execution,
          ::Kokkos::subview(res.triplets.matrix.data(), range),
          host(constraints.values.data(), extra));
        // Every expanded contribution is kept, zero or not, as on the host.
        std::vector<Segment> segments;
        if (res.triplets.size > 0)
          segments.push_back({0, res.triplets.size, 1, 0, true});
        res.triplets.ranks = rank(segments, "Rodin expanded segments");
        return res;
      }

      /**
       * @brief Downloads the number of lifts before the contributions of
       * every matrix.
       * @param starts Number of lifts before every rank.
       * @param offsets Rank of the first contribution of every matrix,
       * followed by their count.
       * @return Number of lifts before the contributions of every entry of
       * @p offsets.
       */
      std::vector<Index> download(
        const Array<Index>& starts, const std::vector<size_t>& offsets) const
      {
        std::vector<Index> res(offsets.size());
        for (size_t k = 0; k < offsets.size(); ++k)
        {
          ::Kokkos::deep_copy(m_execution,
            ::Kokkos::View<Index, ::Kokkos::HostSpace,
              ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(&res[k]),
            ::Kokkos::subview(starts, offsets[k]));
        }
        m_execution.fence("Rodin lifts");
        return res;
      }

      /**
       * @brief Adds listed contributions to the right-hand side, in rank
       * order.
       * @param indices Rows of the contributions.
       * @param values Contributions.
       * @param ranks Ranks of the contributions.
       * @param x Expansions on the masters, or null.
       * @param rhs Device right-hand side, incremented.
       */
      void add(const Array<Index>& indices, const Array<Real>& values, const Ranks& ranks,
        const Expansions* x, const Array<Real>& rhs)
      {
        if (x)
        {
          Array<Index> ei;
          Array<Real> ev;
          expand(m_execution, indices, values, ranks, *x, ei, ev);
          add(ei, ev, Ranks(), nullptr, rhs);
          return;
        }
        const auto sorted = ranks.sorting(m_rows);
        const size_t count = indices.extent(0);
        Listed listed{allocate<Index>(count, "Rodin structural vector keys"),
          allocate<Index>(count, "Rodin structural vector positions")};
        ::Kokkos::parallel_for("Rodin::Assembly::Kokkos::structuralVector",
          Range(m_execution, 0, count),
          VectorKeys{indices, sorted, listed.keys, listed.positions});
        const auto groups = groupContributions(listed, sorted, values, count);
        accumulate(m_execution, listed.positions, groups, values, rhs);
      }

      /**
       * @brief Sums the linear integrators and the preassembled vectors into
       * the right-hand side, in the order they were added.
       * @param rhs Device right-hand side, incremented.
       * @param x Expansions on the masters, or null.
       * @param lift Sums the values lifted by a number of first matrices.
       */
      template <class Lift>
      void sum(const Array<Real>& rhs, const Expansions* x, const Lift& lift)
      {
        for (size_t begin = 0; begin < m_vectors.size();)
        {
          lift(m_vectors[begin].matrices);
          if (!m_vectors[begin].integrator)
          {
            const auto& v = m_vectors[begin++];
            const size_t count =
              v.backend ? v.backend->size() : static_cast<size_t>(v.vec->size());
            const auto values = v.backend
              ? v.backend->getView()
              : upload(v.vec->data(), count, "Rodin preassembled vector");
            if (x)
              add(sequence(m_execution, count, v.row), values, Ranks(), x, rhs);
            else
              add(m_execution, values, v.row, rhs);
            continue;
          }
          size_t end = begin;
          while (end < m_vectors.size() && m_vectors[end].integrator &&
            m_vectors[end].matrices == m_vectors[begin].matrices)
            ++end;
          std::vector<Segment> segments;
          std::vector<size_t> counts;
          size_t total = 0;
          for (size_t k = begin; k < end; ++k)
          {
            auto* integrator = m_vectors[k].integrator;
            const size_t count = integrator->prepare();
            counts.push_back(count);
            if (count > 0)
              segments.push_back({total, count, integrator->getSpacing(), 0, false});
            total += count;
          }
          const auto indices = allocate<Index>(total, "Rodin vector indices");
          const auto values = allocate<Real>(total, "Rodin vector values");
          size_t offset = 0;
          ::Rodin::Kokkos::Internal::clearLaunches();
          for (size_t k = begin; k < end; ++k)
          {
            const auto& v = m_vectors[k];
            const size_t count = counts[k - begin];
            v.integrator->fill(values, offset, v.sign, false, indices);
            ::Rodin::Kokkos::Internal::checkLaunches("A device linear integrator");
            if (v.row != 0)
              shift(m_execution, offset, count, v.row, 0, indices, Array<Index>());
            offset += count;
          }
          add(indices, values, rank(segments, "Rodin vector segments"), x, rhs);
          begin = end;
        }
      }

      ExecutionSpace m_execution;
      size_t m_rows;
      size_t m_cols;
      ConstraintsType m_constraints;
      std::vector<Matrix> m_matrices;
      std::vector<Vector> m_vectors;
      Symbolic* m_symbolic;
      size_t m_group;
  };
}

#endif
