/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_SCATTERMAP_H
#define RODIN_ASSEMBLY_SCATTERMAP_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <type_traits>
#include <vector>

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

#include "Rodin/Types.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"

namespace Rodin::Assembly
{
  /**
   * @brief Cached map from polytope-local matrix entries into the value array
   * of an assembled Math::SparseMatrix.
   *
   * The sparsity pattern of a named form is fixed by the mesh, the DOF maps
   * of the two spaces and the polytopes the form integrates over, so it is
   * computed once, on the first assembly, out of the triplet list that
   * assembly produces anyway. Every subsequent assembly zeroes the value array
   * and scatters the local matrices straight into it through the cached
   * indices, skipping both the triplet list and Eigen's symbolic pass.
   *
   * The pattern is purely topological: it holds every entry the DOF maps
   * connect, including those whose value happens to vanish. Reassembling a
   * form therefore never changes the sparsity of its operator, which is what
   * makes the value array reusable.
   *
   * # Regions
   *
   * The map iterates whatever the iteration object selects. Of its
   * candidates, a polytope is visited when @c filter accepts it and, if an
   * attribute set is given, when its attribute belongs to that set; this is
   * the same selection the triplet assembly of an integrator makes. A
   * candidate that is not visited owns an empty range of the index array, so
   * the offsets stay indexed by candidate and need no remapping.
   *
   * # Validity
   *
   * The cached indices are only valid for the DOF maps and the selection they
   * were built from, and the shape of the operator identifies neither: an
   * edge flip or a renumbering can leave the row count, the column count, the
   * nonzero count and the candidate count all unchanged while moving every
   * entry, and so can a change of attribute set. The cache is therefore keyed
   * on a fingerprint of the visited polytopes and their maps.
   *
   * The fingerprint is accumulated inside the scatter loop rather than in a
   * pass of its own: the loop already reads both DOF arrays, so folding the
   * hash into it costs one mix per degree of freedom instead of a second
   * traversal, which measured 2-3x on the reassembly benchmarks. It is
   * compared once the loop is done, and a mismatch discards the values and
   * rebuilds the pattern from triplets, so nothing assembled against a stale
   * pattern is ever observable. Before a polytope scatters, its local size is
   * checked against the range the cache reserved for it, so that a polytope
   * the cache never visited cannot write past its range before the
   * fingerprint gets a chance to reject the pattern. The per-polytope hashes
   * are combined by addition, which keeps the fingerprint independent of the
   * order the threads happen to visit the polytopes in.
   */
  template <class Scalar>
  class ScatterMap
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Assembled operator type.
      using MatrixType = Math::SparseMatrix<ScalarType>;

      /**
       * @brief Assembles @p out sequentially.
       *
       * @param[in,out] out Operator receiving the assembly.
       * @param[in] prototype Local kernel, copied once for the assembly.
       * @param[in] trialFES Trial space.
       * @param[in] testFES Test space.
       * @param[in] seq Iteration over the region to integrate.
       * @param[in] attributes Attributes to restrict to, or empty for all.
       *
       * @tparam KernelType Local kernel, providing
       * @c compute(Math::Matrix<ScalarType>&, const Geometry::Polytope&).
       */
      template <class KernelType, class IterationType, class TrialFES, class TestFES>
      void assemble(MatrixType& out, const KernelType& prototype,
        const TrialFES& trialFES, const TestFES& testFES, const IterationType& seq,
        const FlatSet<Geometry::Attribute>& attributes) const
      {
        const size_t d = seq.getDimension();
        const Index count = seq.getCount();
        const auto& mesh = trialFES.getMesh();

        if (isValid(out, trialFES, testFES, count))
        {
          auto* const values = out.valuePtr();
          std::fill(values, values + out.nonZeros(), ScalarType(0));
          KernelType kernel(prototype);
          Math::Matrix<ScalarType> local;
          size_t fingerprint = 0;
          bool stale = false;
          for (Index i = 0; i < count; ++i)
          {
            if (!isVisited(seq, mesh, attributes, d, i))
              continue;
            const auto polytope = seq.getPolytope(i);
            kernel.compute(local, polytope);
            const auto& rows = testFES.getDOFs(d, i);
            const auto& cols = trialFES.getDOFs(d, i);
            fingerprint += hashCell(d, i, rows, cols);
            size_t k = m_offsets[i];
            if (m_offsets[i + 1] - k !=
              static_cast<size_t>(rows.size()) * static_cast<size_t>(cols.size()))
            {
              stale = true;
              break;
            }
            for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
              for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
                values[m_indices[k++]] += Math::conj(local(r, c));
          }
          if (!stale && fingerprint == m_fingerprint)
            return;
          // The maps or the selection moved under the cached pattern, so the
          // values just scattered are meaningless. Fall through and rebuild.
        }

        m_offsets.assign(static_cast<size_t>(count) + 1, 0);
        std::vector<Eigen::Triplet<ScalarType>> triplets;
        triplets.reserve(getCapacity(trialFES, testFES));
        KernelType kernel(prototype);
        Math::Matrix<ScalarType> local;
        size_t fingerprint = 0;
        for (Index i = 0; i < count; ++i)
        {
          if (!isVisited(seq, mesh, attributes, d, i))
            continue;
          const auto polytope = seq.getPolytope(i);
          kernel.compute(local, polytope);
          fingerprint += emplace(triplets, local, trialFES, testFES, d, i);
          m_offsets[i + 1] = static_cast<size_t>(local.size());
        }
        setFromTriplets(out, triplets, trialFES, testFES);
        build(out, trialFES, testFES, d, count, fingerprint);
      }

#ifdef RODIN_USE_OPENMP
      /**
       * @brief Assembles @p out with OpenMP.
       *
       * Each thread owns a copy of the kernel and a contiguous range of
       * candidates, and scatters into the shared value array; the additions
       * that collide on a shared degree of freedom are made atomically.
       *
       * @param[in,out] out Operator receiving the assembly.
       * @param[in] prototype Local kernel, copied once per thread.
       * @param[in] trialFES Trial space.
       * @param[in] testFES Test space.
       * @param[in] seq Iteration over the region to integrate.
       * @param[in] attributes Attributes to restrict to, or empty for all.
       * @param[in] threadCount Number of threads.
       *
       * @tparam KernelType Local kernel, providing
       * @c compute(Math::Matrix<ScalarType>&, const Geometry::Polytope&).
       */
      template <class KernelType, class IterationType, class TrialFES, class TestFES>
      void assemble(MatrixType& out, const KernelType& prototype,
        const TrialFES& trialFES, const TestFES& testFES, const IterationType& seq,
        const FlatSet<Geometry::Attribute>& attributes, int threadCount) const
      {
        const size_t d = seq.getDimension();
        const Index count = seq.getCount();
        const auto& mesh = trialFES.getMesh();

        if (isValid(out, trialFES, testFES, count))
        {
          auto* const values = out.valuePtr();
          const Index nonZeroCount = out.nonZeros();
          size_t fingerprint = 0;
          bool stale = false;
#pragma omp parallel num_threads(threadCount) reduction(+ : fingerprint)                 \
  reduction(|| : stale)
          {
#pragma omp for
            for (Index i = 0; i < nonZeroCount; ++i)
              values[i] = ScalarType(0);

            KernelType kernel(prototype);
            Math::Matrix<ScalarType> local;
#pragma omp for
            for (Index i = 0; i < count; ++i)
            {
              if (!isVisited(seq, mesh, attributes, d, i))
                continue;
              const auto polytope = seq.getPolytope(i);
              kernel.compute(local, polytope);
              const auto& rows = testFES.getDOFs(d, i);
              const auto& cols = trialFES.getDOFs(d, i);
              fingerprint += hashCell(d, i, rows, cols);
              size_t k = m_offsets[i];
              if (m_offsets[i + 1] - k !=
                static_cast<size_t>(rows.size()) * static_cast<size_t>(cols.size()))
              {
                stale = true;
                continue;
              }
              for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
              {
                for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
                {
                  const ScalarType s = Math::conj(local(r, c));
                  add(values[m_indices[k++]], s);
                }
              }
            }
          }
          if (!stale && fingerprint == m_fingerprint)
            return;
          // The maps or the selection moved under the cached pattern, so the
          // values just scattered are meaningless. Fall through and rebuild.
        }

        m_offsets.assign(static_cast<size_t>(count) + 1, 0);
        const size_t capacity = getCapacity(trialFES, testFES);
        std::vector<std::vector<Eigen::Triplet<ScalarType>>> chunks(
          static_cast<size_t>(threadCount));
        size_t fingerprint = 0;

#pragma omp parallel num_threads(threadCount) reduction(+ : fingerprint)
        {
          KernelType kernel(prototype);
          Math::Matrix<ScalarType> local;
          auto& triplets = chunks[static_cast<size_t>(omp_get_thread_num())];
          triplets.reserve(capacity / std::max(1, threadCount));
#pragma omp for
          for (Index i = 0; i < count; ++i)
          {
            if (!isVisited(seq, mesh, attributes, d, i))
              continue;
            const auto polytope = seq.getPolytope(i);
            kernel.compute(local, polytope);
            fingerprint += emplace(triplets, local, trialFES, testFES, d, i);
            m_offsets[i + 1] = static_cast<size_t>(local.size());
          }
        }

        std::vector<Eigen::Triplet<ScalarType>> triplets;
        triplets.reserve(capacity);
        for (auto& chunk : chunks)
          triplets.insert(triplets.end(), std::make_move_iterator(chunk.begin()),
            std::make_move_iterator(chunk.end()));
        setFromTriplets(out, triplets, trialFES, testFES);
        build(out, trialFES, testFES, d, count, fingerprint);
      }
#endif

    private:
#ifdef RODIN_USE_OPENMP
      /**
       * @brief Atomically accumulates @p s into @p value.
       *
       * Only the OpenMP scatter needs this, and the @c omp @c atomic
       * directives are not valid outside an OpenMP build.
       */
      static void add(ScalarType& value, const ScalarType& s)
      {
        if constexpr (std::is_floating_point_v<ScalarType>)
        {
#pragma omp atomic
          value += s;
        }
        else
        {
          // std::complex<T> is layout-compatible with T[2].
          auto* const parts = reinterpret_cast<typename ScalarType::value_type*>(&value);
#pragma omp atomic
          parts[0] += s.real();
#pragma omp atomic
          parts[1] += s.imag();
        }
      }
#endif

      /**
       * @brief Tests whether candidate @p i of @p seq is integrated.
       *
       * The iteration decides membership in the region; a non-empty
       * @p attributes further restricts to polytopes carrying one of them.
       */
      template <class IterationType, class MeshType>
      static bool isVisited(const IterationType& seq, const MeshType& mesh,
        const FlatSet<Geometry::Attribute>& attributes, size_t d, Index i)
      {
        if (!seq.filter(i))
          return false;
        if (attributes.empty())
          return true;
        const auto attribute = mesh.getAttribute(d, i);
        return attribute && attributes.count(*attribute);
      }

      template <class TrialFES, class TestFES>
      static size_t getCapacity(const TrialFES& trialFES, const TestFES& testFES)
      {
        return testFES.getSize() * std::log(trialFES.getSize());
      }

      /**
       * @brief Appends the entries of one local matrix to @p triplets.
       * @returns The polytope's contribution to the fingerprint.
       */
      template <class TrialFES, class TestFES>
      static size_t emplace(std::vector<Eigen::Triplet<ScalarType>>& triplets,
        const Math::Matrix<ScalarType>& local, const TrialFES& trialFES,
        const TestFES& testFES, size_t d, Index i)
      {
        const auto& rows = testFES.getDOFs(d, i);
        const auto& cols = trialFES.getDOFs(d, i);
        for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
          for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
            triplets.emplace_back(rows(r), cols(c), Math::conj(local(r, c)));
        return hashCell(d, i, rows, cols);
      }

      template <class TrialFES, class TestFES>
      static void setFromTriplets(MatrixType& out,
        const std::vector<Eigen::Triplet<ScalarType>>& triplets, const TrialFES& trialFES,
        const TestFES& testFES)
      {
        out.resize(testFES.getSize(), trialFES.getSize());
        out.setFromTriplets(triplets.begin(), triplets.end());
        out.makeCompressed();
      }

      /**
       * @brief Combines @p value into @p seed.
       *
       * boost::hash_combine's mixer, kept local so that ScatterMap does not
       * depend on Boost.
       */
      static void mix(size_t& seed, size_t value)
      {
        seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
      }

      /**
       * @brief Hashes the entries one visited polytope contributes.
       *
       * Covers the candidate index, the iteration dimension and both DOF
       * arrays, so any rewiring or renumbering of the maps shows up even when
       * every count stays the same, and so does a change in which candidates
       * are visited. The whole-iteration fingerprint is the sum of these,
       * which makes it independent of visiting order.
       */
      template <class Rows, class Cols>
      static size_t hashCell(size_t d, Index i, const Rows& rows, const Cols& cols)
      {
        size_t seed = 0;
        mix(seed, d);
        mix(seed, static_cast<size_t>(i));
        mix(seed, static_cast<size_t>(rows.size()));
        for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
          mix(seed, static_cast<size_t>(rows(r)));
        mix(seed, static_cast<size_t>(cols.size()));
        for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
          mix(seed, static_cast<size_t>(cols(c)));
        return seed;
      }

      /**
       * @brief Tests whether the cached pattern still fits @p out.
       *
       * The matrix dimensions and compressed sparse index arrays must match
       * the pattern for which the cached scatter indices were built. Whether
       * it fits the DOF maps and the selection is settled by the fingerprint,
       * which the scatter loop accumulates as it goes.
       */
      template <class TrialFES, class TestFES>
      bool isValid(const MatrixType& out, const TrialFES& trialFES,
        const TestFES& testFES, Index count) const
      {
        if (!m_built || !out.isCompressed() ||
          static_cast<size_t>(out.rows()) != testFES.getSize() ||
          static_cast<size_t>(out.cols()) != trialFES.getSize() ||
          static_cast<size_t>(out.nonZeros()) != m_nonZeroCount ||
          m_offsets.size() != static_cast<size_t>(count) + 1 ||
          m_outerIndices.size() != static_cast<size_t>(out.outerSize()) + 1 ||
          m_innerIndices.size() != static_cast<size_t>(out.nonZeros()))
          return false;
        return std::equal(
                 m_outerIndices.begin(), m_outerIndices.end(), out.outerIndexPtr()) &&
          std::equal(m_innerIndices.begin(), m_innerIndices.end(), out.innerIndexPtr());
      }

      /**
       * @brief Locates every visited local entry in the value array of @p out.
       *
       * Expects @c m_offsets to hold, at <tt>i + 1</tt>, the local size of
       * candidate @c i, zero when it was not visited; turns them into offsets
       * and fills the index array.
       */
      template <class TrialFES, class TestFES>
      void build(const MatrixType& out, const TrialFES& trialFES, const TestFES& testFES,
        size_t d, Index count, size_t fingerprint) const
      {
        const auto* const outer = out.outerIndexPtr();
        const auto* const inner = out.innerIndexPtr();
        for (Index i = 0; i < count; ++i)
          m_offsets[i + 1] += m_offsets[i];
        m_indices.resize(m_offsets.back());
        for (Index i = 0; i < count; ++i)
        {
          if (m_offsets[i + 1] == m_offsets[i])
            continue;
          const auto& rows = testFES.getDOFs(d, i);
          const auto& cols = trialFES.getDOFs(d, i);
          size_t k = m_offsets[i];
          for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
          {
            const auto row = static_cast<typename MatrixType::StorageIndex>(rows(r));
            for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
            {
              // Math::SparseMatrix is column major with sorted inner indices.
              const auto column = cols(c);
              const auto* const begin = inner + outer[column];
              const auto* const end = inner + outer[column + 1];
              const auto* const it = std::lower_bound(begin, end, row);
              assert(it != end && *it == row);
              m_indices[k++] = static_cast<Index>(it - inner);
            }
          }
          assert(k == m_offsets[i + 1]);
        }
        m_nonZeroCount = static_cast<size_t>(out.nonZeros());
        m_outerIndices.assign(
          out.outerIndexPtr(), out.outerIndexPtr() + out.outerSize() + 1);
        m_innerIndices.assign(out.innerIndexPtr(), out.innerIndexPtr() + out.nonZeros());
        m_fingerprint = fingerprint;
        m_built = true;
      }

      mutable bool m_built = false;
      mutable size_t m_nonZeroCount = 0;
      mutable std::vector<typename MatrixType::StorageIndex> m_outerIndices;
      mutable std::vector<typename MatrixType::StorageIndex> m_innerIndices;
      mutable size_t m_fingerprint = 0;
      mutable std::vector<size_t> m_offsets;
      mutable std::vector<Index> m_indices;
  };
}

#endif
