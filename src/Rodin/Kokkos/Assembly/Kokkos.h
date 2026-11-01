/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_KOKKOS_H
#define RODIN_KOKKOS_ASSEMBLY_KOKKOS_H

/**
 * @file
 * @brief Device assembly of bilinear forms, linear forms, problems and
 * Dirichlet conditions on resident Kokkos meshes.
 */

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Mult.h"
#include "Rodin/Variational/UnaryMinus.h"

#include "Rodin/Kokkos/ForwardDecls.h"
#include "Rodin/Kokkos/Variational/Integral.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Finds the device implementation of a local integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lbfi Local integrator, possibly negated or scaled by a real.
   * @param sign Factor of the contributions.
   * @return Device integrator and the factor of its contributions.
   */
  template <class Backend>
  std::pair<const ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>*,
    Real>
  getDeviceIntegrator(
    const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Real sign = Real(1))
  {
    using Negated =
      Variational::UnaryMinus<Variational::LocalBilinearFormIntegratorBase<Real>>;
    using Scaled =
      Variational::Mult<Real, Variational::LocalBilinearFormIntegratorBase<Real>>;
    if (const auto* negated = dynamic_cast<const Negated*>(&lbfi))
      return getDeviceIntegrator<Backend>(negated->getOperand(), -sign);
    if (const auto* scaled = dynamic_cast<const Scaled*>(&lbfi))
      return getDeviceIntegrator<Backend>(scaled->getRHS(), scaled->getLHS() * sign);
    const auto* integrator = dynamic_cast<
      const ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>*>(&lbfi);
    if (!integrator)
      Alert::Exception() << "Local integrator has no device implementation."
                         << Alert::Raise;
    return {integrator, sign};
  }

  /**
   * @brief Finds the device implementation of a global integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param gbfi Global integrator, possibly negated.
   * @param sign Factor of the contributions.
   * @return Device integrator and the factor of its contributions.
   */
  template <class Backend>
  std::pair<const ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>*,
    Real>
  getDeviceIntegrator(
    const Variational::GlobalBilinearFormIntegratorBase<Real>& gbfi, Real sign = Real(1))
  {
    using Negated =
      Variational::UnaryMinus<Variational::GlobalBilinearFormIntegratorBase<Real>>;
    if (const auto* negated = dynamic_cast<const Negated*>(&gbfi))
      return getDeviceIntegrator<Backend>(negated->getOperand(), -sign);
    const auto* integrator = dynamic_cast<
      const ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>*>(&gbfi);
    if (!integrator)
      Alert::Exception() << "Global integrator has no device implementation."
                         << Alert::Raise;
    return {integrator, sign};
  }

  /**
   * @brief Finds the device implementation of a linear integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lfi Integrator, possibly negated or scaled by a real.
   * @param sign Factor of the contributions.
   * @return Device integrator and the factor of its contributions.
   */
  template <class Backend>
  std::pair<const ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>*,
    Real>
  getDeviceIntegrator(
    const Variational::LinearFormIntegratorBase<Real>& lfi, Real sign = Real(1))
  {
    using Negated = Variational::UnaryMinus<Variational::LinearFormIntegratorBase<Real>>;
    using Scaled = Variational::Mult<Real, Variational::LinearFormIntegratorBase<Real>>;
    if (const auto* negated = dynamic_cast<const Negated*>(&lfi))
      return getDeviceIntegrator<Backend>(negated->getOperand(), -sign);
    if (const auto* scaled = dynamic_cast<const Scaled*>(&lfi))
      return getDeviceIntegrator<Backend>(scaled->getRHS(), scaled->getLHS() * sign);
    const auto* integrator = dynamic_cast<
      const ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>*>(&lfi);
    if (!integrator)
      Alert::Exception() << "Integrator has no device implementation." << Alert::Raise;
    return {integrator, sign};
  }

  /**
   * @brief Linear system assembled in device memory with value constraints
   * eliminated, then copied to Eigen.
   *
   * The system is bitwise the one of the host assemblies. The nonzero
   * entries of the local integrators, then of the global integrators, in
   * rank order, followed by the entries of the preassembled forms make the
   * list of matrix
   * contributions. Rows of fixed degrees of freedom are dropped and become
   * identity rows holding the fixed value, and the columns of fixed degrees
   * of freedom move to the right-hand side. The contributions to every
   * coefficient are summed in list order, as
   * Eigen::SparseMatrix::setFromTriplets() sums duplicates, and products
   * are rounded before being summed. The entries stay where the
   * integrators store them, their list order following from their
   * position, and only the compressed matrix and the right-hand side are
   * copied to the host.
   *
   * Identification constraints are not supported.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class DeviceSystem
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space of the kernels.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Memory space of the system.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device array.
      template <class T>
      using Array = ::Kokkos::View<T*, MemorySpace>;

      /// @brief Device bilinear integrator.
      using BilinearType =
        ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>;

      /// @brief Device linear integrator.
      using LinearType =
        ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>;

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

        /// @brief Whether the zero entries are kept, as for a preassembled
        /// operator.
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

        /// @brief Finds the segment holding a position.
        KOKKOS_INLINE_FUNCTION size_t find(size_t e) const
        {
          size_t k = 0;
          while (k < segments.extent(0) && e >= segments(k).offset + segments(k).count)
            ++k;
          return k;
        }

        /// @brief Rank of the entry stored at a position.
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

        /// @brief Number of the polytope of the entry stored at a position.
        KOKKOS_INLINE_FUNCTION Index polytope(size_t e) const
        {
          const size_t k = find(e);
          if (k == segments.extent(0))
            return polytopes - 1;
          const Segment& s = segments(k);
          return s.spacing > 1 ? s.first + (e - s.offset) % s.spacing : s.first;
        }

        /// @brief Sort key of the entry stored at a position.
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

        /// @brief Whether the entry stored at a position is listed.
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
        /// @brief Row indices.
        Array<Index> rows;

        /// @brief Column indices.
        Array<Index> cols;

        /// @brief Contributions.
        Array<Real> values;

        /// @brief Ranks of the entries.
        Ranks ranks;

        /// @brief Number of entries.
        size_t size = 0;
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

        /// @brief Adds counts.
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
        /// @brief Start of every group in the entries, followed by the entry
        /// count.
        Array<Index> starts;

        /// @brief Key of every group.
        Array<Index> keys;

        /// @brief Number of groups.
        size_t size() const
        {
          return keys.extent(0);
        }
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
         */
        KOKKOS_INLINE_FUNCTION Real operator()(Index e) const
        {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
          KOKKOS_IF_ON_DEVICE((return -__dmul_rn(values(e), fixedValues(cols(e)));))
#endif
          return -(values(e) * fixedValues(cols(e)));
        }
      };

      /**
       * @brief Creates an empty system.
       * @param execution Execution space instance of the mesh.
       * @param rows Number of rows, the size of the test spaces.
       * @param cols Number of columns, the size of the trial spaces.
       */
      DeviceSystem(const ExecutionSpace& execution, size_t rows, size_t cols)
        : m_execution(execution),
          m_rows(rows),
          m_cols(cols),
          m_fixed(std::max(rows, cols), 0),
          m_values(std::max(rows, cols), Real(0)),
          m_constrained(false)
      {}

      /**
       * @brief Assembles a Dirichlet condition and fixes its degrees of
       * freedom.
       * @param dbc Value Dirichlet condition.
       * @param offset Offset of the constrained space in the system.
       */
      void fix(Variational::DirichletBCBase<Real>& dbc, size_t offset)
      {
        using ValueDOFs = typename Variational::DirichletBCBase<Real>::ValueDOFs;
        dbc.assemble();
        const auto* dofs = std::get_if<ValueDOFs>(&dbc.getDOFs());
        if (!dofs)
        {
          Alert::Exception()
            << "Identification conditions are not supported on Kokkos meshes."
            << Alert::Raise;
        }
        for (const auto& [local, value] : *dofs)
        {
          m_fixed[offset + local] = 1;
          m_values[offset + local] = value;
          m_constrained = true;
        }
      }

      /**
       * @brief Adds the entries of a local integrator to the matrix.
       * @param lbfi Local integrator, possibly negated or scaled by a real.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi,
        Index row, Index col)
      {
        const auto [integrator, sign] = getDeviceIntegrator<Backend>(lbfi);
        m_matrices.push_back({integrator, nullptr, row, col, sign});
      }

      /**
       * @brief Adds the entries of a global integrator to the matrix.
       * @param gbfi Global integrator, possibly negated.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(const Variational::GlobalBilinearFormIntegratorBase<Real>& gbfi,
        Index row, Index col)
      {
        const auto [integrator, sign] = getDeviceIntegrator<Backend>(gbfi);
        m_matrices.push_back({integrator, nullptr, row, col, sign});
      }

      /**
       * @brief Adds the entries of a preassembled operator to the matrix.
       * @param op Operator, alive until finish() returns.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(const Math::SparseMatrix<Real>& op, Index row, Index col)
      {
        m_matrices.push_back({nullptr, &op, row, col, Real(1)});
      }

      /**
       * @brief Adds the entries of a linear integrator to the right-hand side.
       * @param lfi Integrator, possibly negated or scaled by a real.
       * @param row Offset of its test space.
       * @param sign Factor of the contributions.
       */
      void addVector(
        const Variational::LinearFormIntegratorBase<Real>& lfi, Index row, Real sign)
      {
        const auto [integrator, s] = getDeviceIntegrator<Backend>(lfi, sign);
        m_vectors.push_back({integrator, nullptr, row, s});
      }

      /**
       * @brief Adds a preassembled vector to the right-hand side.
       * @param vec Vector, alive until finish() returns.
       * @param row Offset of its test space.
       */
      void addVector(const Math::Vector<Real>& vec, Index row)
      {
        m_vectors.push_back({nullptr, &vec, row, Real(1)});
      }

      /**
       * @brief Integrates, eliminates the constraints and builds the system.
       * @param A Matrix, built when not null.
       * @param b Right-hand side, built when not null.
       */
      void finish(Math::SparseMatrix<Real>* A, Math::Vector<Real>* b) const
      {
        Array<unsigned char> fixed;
        Array<Real> fixedValues;
        if (m_constrained)
        {
          fixed = upload(m_fixed, "Rodin fixed DOFs");
          fixedValues = upload(m_values, "Rodin fixed values");
        }
        // The identity rows come last, ranked by position.
        const size_t diagonal = A && m_constrained ? std::min(m_rows, m_cols) : 0;
        const bool eliminate = b && m_constrained;
        const auto coefficientRanks = [&](const Ranks& ranks) {
          return ranks.sorting(m_rows * m_cols);
        };
        Triplets triplets;
        Listed coefficients, eliminated;
        if (A || eliminate)
        {
          triplets = listMatrix();
          listContributions(m_execution, m_rows, diagonal, triplets, fixed,
            coefficientRanks(triplets.ranks), triplets.ranks.sorting(m_rows),
            A != nullptr, eliminate, coefficients, eliminated);
        }
        if (A)
          compress(triplets, coefficientRanks(triplets.ranks), coefficients, *A);
        if (b)
        {
          const Array<Real> rhs(
            ::Kokkos::view_alloc(m_execution, std::string("Rodin right-hand side")),
            m_rows);
          sumVectors(rhs);
          if (m_constrained)
          {
            accumulate(m_execution, eliminated.positions,
              group(m_execution, eliminated.keys, eliminated.positions,
                triplets.ranks.sorting(m_rows)),
              Eliminated{triplets.values, triplets.cols, fixedValues}, rhs);
            imposeFixed(m_execution, m_rows, fixed, fixedValues, rhs);
          }
          b->resize(m_rows);
          ::Kokkos::deep_copy(m_execution, host(b->data(), m_rows), rhs);
        }
        m_execution.fence("Rodin linear system");
      }

      /**
       * @brief Adds the offsets of their spaces to the indices of entries.
       * @param cols Column indices, empty for a vector.
       *
       * Public as the enclosing function of device lambdas.
       */
      static void shift(const ExecutionSpace& execution, size_t offset, size_t count,
        Index row, Index col, const Array<Index>& rows, const Array<Index>& cols)
      {
        const bool matrix = cols.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::shift",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, offset, offset + count),
          KOKKOS_LAMBDA(const size_t e) {
            rows(e) += row;
            if (matrix)
              cols(e) += col;
          });
      }

      /**
       * @brief Sorts keys with the positions of their entries.
       * @param keys Sort keys, sorted on return.
       * @param positions Positions of the entries, permuted as the keys on
       * return, the positions of equal keys being in rank order.
       * @param ranks Ranks of the entries, with the layout of the sort keys.
       * @returns Groups of equal keys.
       *
       * The order among equal keys does not depend on the sort being stable.
       *
       * Public as the enclosing function of device lambdas.
       */
      static Groups group(const ExecutionSpace& execution, const Array<Index>& keys,
        const Array<Index>& positions, const Ranks& ranks)
      {
        const size_t count = keys.extent(0);
        if (count > 0)
          ::Kokkos::Experimental::sort_by_key(execution, keys, positions);
        const unsigned shift = ranks.shift;
        const auto head = KOKKOS_LAMBDA(const size_t i)
        {
          return i == 0 || (keys(i) >> shift) != (keys(i - 1) >> shift);
        };
        const auto policy = ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count);
        size_t distinct = 0;
        size_t unordered = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::groups", policy,
          KOKKOS_LAMBDA(const size_t i, size_t& d, size_t& u) {
            if (head(i))
              ++d;
            else
              u += ranks(positions(i)) < ranks(positions(i - 1));
          },
          ::Kokkos::Sum<size_t>(distinct), ::Kokkos::Sum<size_t>(unordered));
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
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::starts", policy,
          KOKKOS_LAMBDA(const size_t i, size_t& update, const bool final) {
            if (head(i))
            {
              if (final)
              {
                starts(update) = i;
                groupKeys(update) = keys(i) >> shift;
              }
              ++update;
            }
          });
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(starts, distinct), count);
        if (unordered == 0)
          return res;

        // Groups are small, a few polytopes sharing a pair of degrees of freedom.
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::order",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, distinct),
          KOKKOS_LAMBDA(const size_t u) {
            for (Index i = starts(u) + 1; i < starts(u + 1); ++i)
            {
              const Index p = positions(i);
              const Index r = ranks(p);
              Index j = i;
              for (; j > starts(u) && ranks(positions(j - 1)) > r; --j)
                positions(j) = positions(j - 1);
              positions(j) = p;
            }
          });
        return res;
      }

      /**
       * @brief Lists the matrix contributions kept on their coefficients
       * with the identity rows, and the contributions moved to the
       * right-hand side by the columns of the fixed degrees of freedom.
       * @param diagonal Number of identity rows listed.
       * @param coefficientRanks Ranks of the contributions, with the layout
       * of the coefficient sort keys.
       * @param eliminatedRanks Ranks of the contributions, with the layout
       * of the row sort keys.
       * @param listCoefficients Whether the coefficients are listed.
       * @param listEliminated Whether the eliminated contributions are listed.
       * @param coefficients Receives the coefficients, as column major sort
       * keys, the diagonal of the fixed row i being at position size + i.
       * @param eliminated Receives the eliminated contributions, as row sort
       * keys.
       *
       * Public as the enclosing function of device lambdas.
       */
      static void listContributions(const ExecutionSpace& execution, size_t rowCount,
        size_t diagonal, const Triplets& triplets, const Array<unsigned char>& fixed,
        const Ranks& coefficientRanks, const Ranks& eliminatedRanks,
        bool listCoefficients, bool listEliminated, Listed& coefficients,
        Listed& eliminated)
      {
        const size_t size = triplets.size;
        const auto rows = triplets.rows;
        const auto cols = triplets.cols;
        const auto values = triplets.values;
        const bool constrained = fixed.extent(0) > 0;
        const auto classify = KOKKOS_LAMBDA(const size_t e)->Counts
        {
          Counts c;
          if (e >= size)
          {
            c.coefficients = fixed(e - size) != 0;
            return c;
          }
          if (!coefficientRanks.kept(values, e))
            return c;
          const bool row = constrained && fixed(rows(e));
          const bool col = constrained && fixed(cols(e));
          c.coefficients = listCoefficients && !row && !col;
          c.eliminated = listEliminated && !row && col;
          return c;
        };
        const auto policy =
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, size + diagonal);
        Counts count;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::countContributions", policy,
          KOKKOS_LAMBDA(const size_t e, Counts& sum) { sum += classify(e); }, count);
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
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::listContributions", policy,
          KOKKOS_LAMBDA(const size_t e, Counts& update, const bool final) {
            const Counts c = classify(e);
            if (final && c.coefficients)
            {
              const Index row = e < size ? rows(e) : e - size;
              const Index col = e < size ? cols(e) : e - size;
              ck(update.coefficients) = coefficientRanks.sortKey(col * rowCount + row, e);
              cp(update.coefficients) = e;
            }
            if (final && c.eliminated)
            {
              ek(update.eliminated) = eliminatedRanks.sortKey(rows(e), e);
              ep(update.eliminated) = e;
            }
            update += c;
          });
      }

      /**
       * @brief Sums the contributions of every coefficient into the
       * compressed column storage.
       *
       * Public as the enclosing function of device lambdas.
       */
      static void sumCoefficients(const ExecutionSpace& execution, size_t rowCount,
        size_t colCount, const Triplets& triplets, const Array<Index>& positions,
        const Groups& groups, const Array<Math::SparseIndex>& outer,
        const Array<Math::SparseIndex>& inner, const Array<Real>& data)
      {
        const size_t size = triplets.size;
        const auto values = triplets.values;
        const auto starts = groups.starts;
        const auto keys = groups.keys;
        const size_t nnz = data.extent(0);
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::sumCoefficients",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, nnz),
          KOKKOS_LAMBDA(const size_t u) {
            const auto value = [&](const Index e) {
              return e < size ? values(e) : Real(1);
            };
            Real sum = value(positions(starts(u)));
            for (Index i = starts(u) + 1; i < starts(u + 1); ++i)
              sum += value(positions(i));
            data(u) = sum;
            inner(u) = static_cast<Math::SparseIndex>(keys(u) % rowCount);
          });
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::columns",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, nnz + 1),
          KOKKOS_LAMBDA(const size_t u) {
            const Index col = u < nnz ? keys(u) / rowCount : colCount;
            const Index begin = u > 0 ? keys(u - 1) / rowCount + 1 : 0;
            for (Index j = begin; j <= col; ++j)
              outer(j) = static_cast<Math::SparseIndex>(u);
          });
      }

      /**
       * @brief Lists the nonzero entries of a vector.
       * @param ranks Ranks of the entries, with the layout of the sort keys.
       * @param keys Receives the sort keys of their indices.
       * @param positions Receives their positions.
       *
       * Public as the enclosing function of device lambdas.
       */
      static void listEntries(const ExecutionSpace& execution,
        const Array<Index>& indices, const Array<Real>& values, const Ranks& ranks,
        Array<Index>& keys, Array<Index>& positions)
      {
        const auto policy =
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, indices.extent(0));
        size_t count = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::countEntries", policy,
          KOKKOS_LAMBDA(const size_t e, size_t& sum) { sum += values(e) != Real(0); },
          count);
        keys = Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                              std::string("Rodin vector keys")),
          count);
        positions =
          Array<Index>(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                         std::string("Rodin vector positions")),
            count);
        const auto k = keys;
        const auto p = positions;
        ::Kokkos::parallel_scan(
          "Rodin::Assembly::Kokkos::listEntries", policy,
          KOKKOS_LAMBDA(const size_t e, size_t& update, const bool final) {
            if (values(e) != Real(0))
            {
              if (final)
              {
                k(update) = ranks.sortKey(indices(e), e);
                p(update) = e;
              }
              ++update;
            }
          });
      }

      /**
       * @brief Adds grouped contributions to the right-hand side, in rank
       * order.
       * @param contributions Contribution of an entry, by position.
       *
       * Public as the enclosing function of device lambdas.
       */
      template <class Contribution>
      static void accumulate(const ExecutionSpace& execution,
        const Array<Index>& positions, const Groups& groups,
        const Contribution& contributions, const Array<Real>& rhs)
      {
        const auto starts = groups.starts;
        const auto keys = groups.keys;
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::accumulate",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, groups.size()),
          KOKKOS_LAMBDA(const size_t u) {
            const Index i = keys(u);
            Real x = rhs(i);
            for (Index k = starts(u); k < starts(u + 1); ++k)
              x += contributions(positions(k));
            rhs(i) = x;
          });
      }

      /**
       * @brief Adds a preassembled vector to the right-hand side.
       *
       * Public as the enclosing function of device lambdas.
       */
      static void addDense(const ExecutionSpace& execution, const Array<Real>& values,
        Index row, const Array<Real>& rhs)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::addDense",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, values.extent(0)),
          KOKKOS_LAMBDA(const size_t i) { rhs(row + i) += values(i); });
      }

      /**
       * @brief Writes the fixed values on the fixed rows.
       *
       * Public as the enclosing function of device lambdas.
       */
      static void imposeFixed(const ExecutionSpace& execution, size_t count,
        const Array<unsigned char>& fixed, const Array<Real>& values,
        const Array<Real>& rhs)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::imposeFixed",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          KOKKOS_LAMBDA(const size_t i) {
            if (fixed(i))
              rhs(i) = values(i);
          });
      }

    private:
      struct Matrix
      {
        const BilinearType* integrator;
        const Math::SparseMatrix<Real>* op;
        Index row;
        Index col;
        Real sign;
      };

      struct Vector
      {
        const LinearType* integrator;
        const Math::Vector<Real>* vec;
        Index row;
        Real sign;
      };

      template <class T>
      static ::Kokkos::View<T*, ::Kokkos::HostSpace,
        ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>
      host(T* data, size_t size)
      {
        return ::Kokkos::View<T*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(data, size);
      }

      template <class T>
      Array<T> allocate(size_t size, const char* label) const
      {
        return Array<T>(::Kokkos::view_alloc(
                          m_execution, ::Kokkos::WithoutInitializing, std::string(label)),
          size);
      }

      template <class T>
      Array<T> upload(const T* data, size_t size, const char* label) const
      {
        const auto device = allocate<T>(size, label);
        ::Kokkos::deep_copy(m_execution, device, host(data, size));
        // The host buffer may be released before the copy.
        m_execution.fence("Rodin upload");
        return device;
      }

      template <class T>
      Array<T> upload(const std::vector<T>& data, const char* label) const
      {
        return upload(data.data(), data.size(), label);
      }

      /**
       * @brief Integrates the local integrators and lists the matrix
       * contributions.
       */
      Triplets listMatrix() const
      {
        std::vector<Segment> segments;
        std::vector<size_t> counts;
        size_t total = 0;
        for (const auto& m : m_matrices)
        {
          const size_t count = m.integrator ? m.integrator->prepare() : m.op->nonZeros();
          counts.push_back(count);
          if (count > 0)
          {
            segments.push_back({total, count,
              m.integrator ? m.integrator->getSpacing() : size_t(1), 0, !m.integrator});
          }
          total += count;
        }
        Triplets triplets;
        triplets.size = total;
        triplets.rows = allocate<Index>(total, "Rodin matrix rows");
        triplets.cols = allocate<Index>(total, "Rodin matrix columns");
        triplets.values = allocate<Real>(total, "Rodin matrix values");
        size_t offset = 0;
        for (size_t k = 0; k < m_matrices.size(); ++k)
        {
          const auto& m = m_matrices[k];
          const size_t count = counts[k];
          if (m.integrator)
          {
            m.integrator->fill(
              triplets.values, offset, m.sign, false, triplets.rows, triplets.cols);
            if (m.row != 0 || m.col != 0)
              shift(
                m_execution, offset, count, m.row, m.col, triplets.rows, triplets.cols);
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
            for (typename Math::SparseMatrix<Real>::InnerIterator it(*m.op, j); it; ++it)
            {
              rows.push_back(m.row + it.row());
              cols.push_back(m.col + it.col());
              values.push_back(it.value());
            }
          }
          const auto range = std::make_pair(offset, offset + count);
          ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(triplets.rows, range),
            host(rows.data(), rows.size()));
          ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(triplets.cols, range),
            host(cols.data(), cols.size()));
          ::Kokkos::deep_copy(m_execution, ::Kokkos::subview(triplets.values, range),
            host(values.data(), values.size()));
          m_execution.fence("Rodin preassembled operator");
          offset += count;
        }
        triplets.ranks = rank(segments, "Rodin matrix segments");
        return triplets;
      }

      /// @brief Numbers the polytopes of segments and uploads them.
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
       * @param ranks Ranks of the contributions, with the layout of the
       * coefficient sort keys.
       */
      void compress(const Triplets& triplets, const Ranks& ranks, const Listed& listed,
        Math::SparseMatrix<Real>& A) const
      {
        A.resize(m_rows, m_cols);
        if (m_rows == 0 || m_cols == 0)
          return;
        const auto groups = group(m_execution, listed.keys, listed.positions, ranks);
        const size_t nnz = groups.size();
        if (nnz > static_cast<size_t>(std::numeric_limits<Math::SparseIndex>::max()))
        {
          Alert::Exception() << "The matrix is too large for the sparse index type."
                             << Alert::Raise;
        }
        const auto outer = allocate<Math::SparseIndex>(m_cols + 1, "Rodin matrix outer");
        const auto inner = allocate<Math::SparseIndex>(nnz, "Rodin matrix inner");
        const auto data = allocate<Real>(nnz, "Rodin matrix data");
        sumCoefficients(m_execution, m_rows, m_cols, triplets, listed.positions, groups,
          outer, inner, data);
        A.resizeNonZeros(static_cast<Eigen::Index>(nnz));
        ::Kokkos::deep_copy(m_execution, host(A.outerIndexPtr(), m_cols + 1), outer);
        ::Kokkos::deep_copy(m_execution, host(A.innerIndexPtr(), nnz), inner);
        ::Kokkos::deep_copy(m_execution, host(A.valuePtr(), nnz), data);
      }

      /**
       * @brief Sums the linear integrators and the preassembled vectors into
       * the right-hand side, in the order they were added.
       */
      void sumVectors(const Array<Real>& rhs) const
      {
        for (size_t begin = 0; begin < m_vectors.size();)
        {
          if (!m_vectors[begin].integrator)
          {
            const auto& v = m_vectors[begin++];
            addDense(m_execution,
              upload(v.vec->data(), static_cast<size_t>(v.vec->size()),
                "Rodin preassembled vector"),
              v.row, rhs);
            continue;
          }
          size_t end = begin;
          while (end < m_vectors.size() && m_vectors[end].integrator)
            ++end;
          std::vector<Segment> segments;
          std::vector<size_t> counts;
          size_t total = 0;
          for (size_t k = begin; k < end; ++k)
          {
            const auto* integrator = m_vectors[k].integrator;
            const size_t count = integrator->prepare();
            counts.push_back(count);
            if (count > 0)
              segments.push_back({total, count, integrator->getSpacing(), 0, false});
            total += count;
          }
          const auto indices = allocate<Index>(total, "Rodin vector indices");
          const auto values = allocate<Real>(total, "Rodin vector values");
          size_t offset = 0;
          for (size_t k = begin; k < end; ++k)
          {
            const auto& v = m_vectors[k];
            const size_t count = counts[k - begin];
            v.integrator->fill(values, offset, v.sign, false, indices);
            if (v.row != 0)
              shift(m_execution, offset, count, v.row, 0, indices, Array<Index>());
            offset += count;
          }
          const auto ranks = rank(segments, "Rodin vector segments").sorting(m_rows);
          Array<Index> keys, positions;
          listEntries(m_execution, indices, values, ranks, keys, positions);
          accumulate(m_execution, positions, group(m_execution, keys, positions, ranks),
            values, rhs);
          begin = end;
        }
      }

      ExecutionSpace m_execution;
      size_t m_rows;
      size_t m_cols;
      std::vector<unsigned char> m_fixed;
      std::vector<Real> m_values;
      bool m_constrained;
      std::vector<Matrix> m_matrices;
      std::vector<Vector> m_vectors;
  };

  /**
   * @brief Rejects the problem terms without a device implementation.
   * @param pb Problem body.
   */
  template <class ProblemBody>
  void requireDeviceProblem(const ProblemBody& pb)
  {
    if (pb.getPBCs().size() > 0)
    {
      Alert::Exception() << "Periodic conditions are not supported on Kokkos meshes."
                         << Alert::Raise;
    }
  }
}

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a sparse bilinear form on a Kokkos mesh.
   *
   * Every local integrator of the form, possibly negated, must implement
   * Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator. Each one
   * integrates its cells on the device, and the matrix is summed and
   * compressed on the device before being copied to the host.
   *
   * Global integrators are integrated after the local ones, as on the host.
   *
   * @tparam Solution Solution variable type.
   * @tparam TrialFES Trial space on a resident Kokkos mesh.
   * @tparam TestFES Test space on the same resident Kokkos mesh.
   */
  template <class Solution, class TrialFES, class TestFES>
  class Kokkos<Math::SparseMatrix<Real>,
    Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>>
    final
    : public AssemblyBase<Math::SparseMatrix<Real>,
        Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>>
  {
    public:
      /// @brief Operator type.
      using OperatorType = Math::SparseMatrix<Real>;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Integrates the form on the device and builds its matrix.
       * @param res Storage for the assembled operator.
       * @param input Input data.
       */
      void execute(OperatorType& res, const InputType& input) override
      {
        Internal::DeviceSystem<BackendType> system(
          input.getTrialFES().getMesh().getContext().getExecutionSpace(),
          input.getTestFES().getSize(), input.getTrialFES().getSize());
        for (const auto& lbfi : input.getLocalBFIs())
          system.addMatrix(lbfi, 0, 0);
        for (const auto& gbfi : input.getGlobalBFIs())
          system.addMatrix(gbfi, 0, 0);
        system.finish(&res, nullptr);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a linear form on a Kokkos mesh.
   *
   * Every integrator of the form, possibly negated, must implement
   * Rodin::Kokkos::Variational::DeviceLinearFormIntegrator. Each one
   * integrates its cells on the device, and the vector is summed on the
   * device before being copied to the host.
   *
   * @tparam FES Test space on a resident Kokkos mesh.
   */
  template <class FES>
  class Kokkos<Math::Vector<Real>, Variational::LinearForm<FES, Math::Vector<Real>>> final
    : public AssemblyBase<Math::Vector<Real>,
        Variational::LinearForm<FES, Math::Vector<Real>>>
  {
    public:
      /// @brief Vector type.
      using VectorType = Math::Vector<Real>;

      /// @brief Linear form type.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<VectorType, LinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the test mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<FES>::MeshType>::ContextType::BackendType;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Assigns the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other) = default;

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the form on the device and builds its vector.
       * @param res Storage for the assembled vector.
       * @param input Input data.
       */
      void execute(VectorType& res, const InputType& input) override
      {
        Internal::DeviceSystem<BackendType> system(
          input.getFES().getMesh().getContext().getExecutionSpace(),
          input.getFES().getSize(), 0);
        for (const auto& lfi : input.getLFIs())
          system.addVector(lfi, 0, Real(1));
        system.finish(nullptr, &res);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a single field problem on a Kokkos mesh.
   *
   * The local integrators and linear integrators of the problem, possibly
   * negated, are integrated on the device. Their entries, the preassembled
   * forms and the value Dirichlet conditions are combined on the device,
   * eliminating the constrained degrees of freedom as the host assemblies
   * do, and only the system is copied to the host.
   *
   * Periodic conditions and identification conditions are not supported.
   *
   * @tparam TrialFunction Trial function on a resident Kokkos mesh.
   * @tparam TestFunction Test function on the same mesh.
   */
  template <class TrialFunction, class TestFunction>
  class Kokkos<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
    Variational::Problem<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
      TrialFunction, TestFunction>>
    final : public AssemblyBase<
              Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
              Variational::Problem<
                Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
                TrialFunction, TestFunction>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType =
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

      /// @brief Problem type.
      using ProblemType =
        Variational::Problem<LinearSystemType, TrialFunction, TestFunction>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the trial mesh context.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<TrialFunction>::FESType>::MeshType>::ContextType::
        BackendType;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Assigns the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other) = default;

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
        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        Internal::requireDeviceProblem(pb);
        Internal::DeviceSystem<BackendType> system(
          u.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(),
          v.getFiniteElementSpace().getSize(), u.getFiniteElementSpace().getSize());
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() == u.getUUID())
            system.fix(dbc, 0);
        }
        for (const auto& lbfi : pb.getLocalBFIs())
          system.addMatrix(lbfi, 0, 0);
        for (const auto& gbfi : pb.getGlobalBFIs())
          system.addMatrix(gbfi, 0, 0);
        for (const auto& bf : pb.getBFs())
          system.addMatrix(bf.getOperator(), 0, 0);
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        for (const auto& lfi : pb.getLFIs())
          system.addVector(lfi, 0, Real(-1));
        for (const auto& lf : pb.getLFs())
          system.addVector(lf.getVector(), 0);
        system.finish(&axb.getOperator(), &axb.getVector());
      }

      /**
       * @brief Assembles the full system and keeps only one side.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        LinearSystemType scratch;
        execute(scratch, input);
        if (target == Variational::AssemblyTarget::LHS)
          axb.getOperator() = std::move(scratch.getOperator());
        else
          axb.getVector() = std::move(scratch.getVector());
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a mixed problem on a Kokkos mesh.
   *
   * Each field occupies a block of the system, at the offsets given by the
   * problem input. Integrators are integrated on the device and their entries
   * are shifted to the blocks of their trial and test functions, then the
   * constraints are eliminated as in the single field assembly.
   *
   * Periodic conditions and identification conditions are not supported.
   */
  template <class U1, class U2, class U3, class... Us>
  class Kokkos<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
    Variational::Problem<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
      U1, U2, U3, Us...>>
    final : public AssemblyBase<
              Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
              Variational::Problem<
                Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>, U1, U2,
                U3, Us...>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType =
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the mesh context of the first field.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType::BackendType;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Assigns the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other) = default;

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
        auto& pb = input.getProblemBody();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets = input.getTestOffsets();
        const auto& trialUUIDMap = input.getTrialUUIDMap();
        const auto& testUUIDMap = input.getTestUUIDMap();
        const auto trialOffset = [&](const auto& uuid) -> Index {
          const auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return static_cast<Index>(trialOffsets[it->second]);
        };
        const auto testOffset = [&](const auto& uuid) -> Index {
          const auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return static_cast<Index>(testOffsets[it->second]);
        };

        Internal::requireDeviceProblem(pb);
        const auto& u1 = input.getTrialFunctions().template get<0>().get();
        Internal::DeviceSystem<BackendType> system(
          u1.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(),
          input.getTotalTestSize(), input.getTotalTrialSize());
        for (auto& dbc : pb.getDBCs())
          system.fix(dbc, trialOffset(dbc.getOperand().getUUID()));
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
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        for (const auto& lfi : pb.getLFIs())
          system.addVector(lfi, testOffset(lfi.getTestFunction().getUUID()), Real(-1));
        for (const auto& lf : pb.getLFs())
          system.addVector(lf.getVector(), testOffset(lf.getTestFunction().getUUID()));
        system.finish(&axb.getOperator(), &axb.getVector());
      }

      /**
       * @brief Assembles the full system and keeps only one side.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        LinearSystemType scratch;
        execute(scratch, input);
        if (target == Variational::AssemblyTarget::LHS)
          axb.getOperator() = std::move(scratch.getOperator());
        else
          axb.getVector() = std::move(scratch.getVector());
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a value Dirichlet condition on a Lagrange
   * Kokkos space.
   *
   * The constrained faces are the faces of codimension one carrying one of
   * the essential attributes or, when no attribute is given, the boundary
   * faces. Each constrained degree of freedom is assigned to its first face,
   * where the value is evaluated at its node. Grid function values are
   * therefore interpolated from the degrees of freedom of that face, which
   * matches the host for continuous data.
   *
   * @tparam Solution Solution variable type.
   * @tparam FES Device space, whose degrees of freedom provide their nodes.
   * @tparam ValueDerived Value function type.
   */
  template <class Solution, class FES, class ValueDerived>
    requires ::Rodin::Kokkos::Variational::DeviceSpace<FES>
  class Kokkos<IndexMap<Real>,
    Variational::DirichletBC<Variational::TrialFunction<Solution, FES>,
      Variational::FunctionBase<ValueDerived>>>
    final : public AssemblyBase<IndexMap<Real>,
              Variational::DirichletBC<Variational::TrialFunction<Solution, FES>,
                Variational::FunctionBase<ValueDerived>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Device degrees of freedom.
      using DOFsType = typename FESType::DOFsType;

      /// @brief Dirichlet condition type.
      using DirichletBCType =
        Variational::DirichletBC<Variational::TrialFunction<Solution, FESType>,
          Variational::FunctionBase<ValueDerived>>;

      /// @brief Parent class.
      using Parent = AssemblyBase<IndexMap<Real>, DirichletBCType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Resident mesh type.
      using MeshType = typename FormLanguage::Traits<FESType>::MeshType;

      /// @brief Execution space of the mesh context.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Memory space of the mesh context.
      using MemorySpace = typename MeshType::MemorySpace;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Evaluates the condition on the device.
       * @param res Constrained degrees of freedom and their values.
       * @param input Boundary condition input.
       */
      void execute(IndexMap<Real>& res, const InputType& input) override
      {
        const auto& fes = input.getOperand().getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& essBdr = input.getEssentialBoundary();
        res.clear();
        const Index d = mesh.getDimension();
        if (d == 0)
          return;
        if (essBdr.empty() && !mesh.getConnectivity().hasIncidence(d - 1, d))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Compute face-to-cell incidence before imposing a boundary condition."
            << Alert::Raise;
        }
        const auto& execution = mesh.getContext().getExecutionSpace();
        const DOFsType dofs = fes.getDeviceDOFs();
        const Index stride = dofs.getStride();
        const Index vdim = fes.getVectorDimension();
        const auto attributes =
          ::Rodin::Kokkos::Variational::Internal::uploadAttributes(execution, essBdr);
        ::Kokkos::View<Index*, MemorySpace> owners(
          ::Kokkos::view_alloc(execution, std::string("Rodin boundary owners")), stride);
        ::Kokkos::deep_copy(execution, owners, std::numeric_limits<Index>::max());
        ::Kokkos::View<Real*, MemorySpace> values(
          ::Kokkos::view_alloc(execution, std::string("Rodin boundary values")),
          stride * vdim);
        const ::Rodin::Kokkos::Variational::DeviceView<
          Variational::FunctionBase<ValueDerived>>
          value(input.getValue());
        evaluate(execution, mesh.getDeviceView(), attributes, owners, values, value, dofs,
          mesh.getFaceCount(), vdim);

        const auto host = ::Kokkos::create_mirror_view(owners);
        const auto hostValues = ::Kokkos::create_mirror_view(values);
        ::Kokkos::deep_copy(execution, host, owners);
        ::Kokkos::deep_copy(execution, hostValues, values);
        execution.fence("Rodin boundary values");
        for (Index r = 0; r < vdim; ++r)
        {
          for (Index i = 0; i < stride; ++i)
          {
            if (host(i) != std::numeric_limits<Index>::max())
              res.insert(res.end(), {i + r * stride, hostValues(i + r * stride)});
          }
        }
      }

      /**
       * @brief Selects the constrained faces and evaluates the value.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param attributes Essential attributes; empty selects the boundary.
       * @param owners First constrained face of every scalar degree of freedom.
       * @param values Value of component r of the degree of freedom i at
       * i + r * stride.
       * @param value Device view of the value function.
       * @param dofs Device degrees of freedom.
       * @param faceCount Number of faces.
       * @param vdim Vector dimension of the space.
       *
       * The counts are given by the host, since the topology of the device
       * mesh is not readable on the host. Public as the enclosing function of
       * device lambdas.
       */
      template <class AttributeView, class ValueView>
      static void evaluate(const ExecutionSpace& execution,
        const ::Rodin::Kokkos::Geometry::MeshView& view, const AttributeView& attributes,
        const ::Kokkos::View<Index*, MemorySpace>& owners,
        const ::Kokkos::View<Real*, MemorySpace>& values, const ValueView& value,
        const DOFsType& dofs, Index faceCount, Index vdim)
      {
        const Index stride = dofs.getStride();
        const auto policy =
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
            execution, 0, faceCount);
        const Index faceDim = view.getDimension() - 1;
        ::Kokkos::parallel_for(
          "Rodin boundary owners", policy, KOKKOS_LAMBDA(Index f) {
            if (!isConstrained(view, attributes, f))
              return;
            const auto indices = dofs.bind(view.getPolytope(faceDim, f));
            for (Index q = 0; q < indices.size(); ++q)
              ::Kokkos::atomic_min(&owners(indices[q]), f);
          });
        ::Kokkos::parallel_for(
          "Rodin boundary values", policy, KOKKOS_LAMBDA(Index f) {
            if (!isConstrained(view, attributes, f))
              return;
            const auto face = view.getPolytope(faceDim, f);
            const auto indices = dofs.bind(face);
            for (Index q = 0; q < indices.size(); ++q)
            {
              const Index dof = indices[q];
              if (owners(dof) != f)
                continue;
              const auto p = face.getPoint(dofs.getNode(face, q));
              store(values, dof, stride, vdim, value.getValue(p));
            }
          });
      }

      /**
       * @brief Stores a scalar value of a degree of freedom.
       * @param values Boundary values.
       * @param dof Scalar degree of freedom.
       * @param stride Offset between components.
       * @param vdim Vector dimension of the space.
       * @param v Value at the node.
       */
      KOKKOS_INLINE_FUNCTION static void store(
        const ::Kokkos::View<Real*, MemorySpace>& values, Index dof,
        [[maybe_unused]] Index stride, [[maybe_unused]] Index vdim, Real v)
      {
        values(dof) = v;
      }

      /**
       * @brief Stores the components of a vector value of a degree of freedom.
       * @param values Boundary values, component r of the degree of freedom i
       * at i + r * stride.
       * @param dof Scalar degree of freedom.
       * @param stride Offset between components.
       * @param vdim Vector dimension of the space.
       * @param v Value at the node.
       */
      KOKKOS_INLINE_FUNCTION static void store(
        const ::Kokkos::View<Real*, MemorySpace>& values, Index dof, Index stride,
        Index vdim, const ::Rodin::Kokkos::Geometry::Vector& v)
      {
        for (Index r = 0; r < vdim; ++r)
          values(dof + r * stride) = v(r);
      }

      /**
       * @brief Tests whether a face is constrained.
       * @param view Device mesh.
       * @param attributes Essential attributes; empty selects the boundary.
       * @param f Face index.
       * @return Whether the face carries the condition.
       */
      template <class AttributeView>
      KOKKOS_INLINE_FUNCTION static bool isConstrained(
        const ::Rodin::Kokkos::Geometry::MeshView& view, const AttributeView& attributes,
        Index f)
      {
        const Index d = view.getDimension();
        if (attributes.extent(0) == 0)
          return view.getConnectivity().getIncidence(d - 1, d, f).size() == 1;
        const auto face = view.getPolytope(d - 1, f);
        return ::Rodin::Kokkos::Variational::isSelected(face, attributes);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @brief Default assembly for a test space on a Kokkos mesh.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class Default<Context::Kokkos<Context::Local, Backend>>
  {
    public:
      /**
       * @brief Device assembly of the object.
       * @tparam LinearAlgebraType Target linear algebra type.
       * @tparam Object Variational form type to assemble.
       */
      template <class LinearAlgebraType, class Object>
      using Type = Kokkos<LinearAlgebraType, Object>;
  };

  /**
   * @brief Default assembly for trial and test spaces on Kokkos meshes.
   * @tparam Backend Kokkos backend of both mesh contexts.
   */
  template <class Backend>
  class Default<Context::Kokkos<Context::Local, Backend>,
    Context::Kokkos<Context::Local, Backend>>
  {
    public:
      /**
       * @brief Device assembly of the object.
       * @tparam LinearAlgebraType Target linear algebra type.
       * @tparam Object Variational form type to assemble.
       */
      template <class LinearAlgebraType, class Object>
      using Type = Kokkos<LinearAlgebraType, Object>;
  };
}

#endif
