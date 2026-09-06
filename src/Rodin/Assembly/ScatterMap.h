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

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"

namespace Rodin::Assembly
{
  /**
   * @brief Cached map from cell-local matrix entries into the value array of
   * an assembled Math::SparseMatrix.
   *
   * The sparsity pattern of a named form is fixed by the mesh and the DOF
   * maps of the two spaces, so it is computed once, on the first assembly,
   * out of the triplet list that assembly produces anyway. Every subsequent
   * assembly zeroes the value array and scatters the cell matrices straight
   * into it through the cached indices, skipping both the triplet list and
   * Eigen's symbolic pass.
   *
   * The pattern is purely topological: it holds every entry the DOF maps
   * connect, including those whose value happens to vanish. Reassembling a
   * form therefore never changes the sparsity of its operator, which is what
   * makes the value array reusable.
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
       * @tparam KernelType Cell kernel, constructible from the two spaces.
       */
      template <class KernelType, class IterationType, class TrialFES, class TestFES>
      void assemble(MatrixType& out, const TrialFES& trialFES, const TestFES& testFES,
        const IterationType& seq, size_t d, Index count) const
      {
        if (!isValid(out, trialFES, testFES, count))
        {
          std::vector<Eigen::Triplet<ScalarType>> triplets;
          triplets.reserve(getCapacity(trialFES, testFES));
          KernelType kernel(trialFES, testFES);
          Math::Matrix<ScalarType> local;
          for (Index i = 0; i < count; ++i)
          {
            const auto polytope = seq.getPolytope(i);
            kernel.compute(local, polytope);
            emplace(triplets, local, trialFES, testFES, d, i);
          }
          setFromTriplets(out, triplets, trialFES, testFES);
          build(out, trialFES, testFES, d, count);
          return;
        }

        auto* const values = out.valuePtr();
        std::fill(values, values + out.nonZeros(), ScalarType(0));
        KernelType kernel(trialFES, testFES);
        Math::Matrix<ScalarType> local;
        for (Index i = 0; i < count; ++i)
        {
          const auto polytope = seq.getPolytope(i);
          kernel.compute(local, polytope);
          const auto& rows = testFES.getDOFs(d, i);
          const auto& cols = trialFES.getDOFs(d, i);
          size_t k = m_offsets[i];
          for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
            for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
              values[m_indices[k++]] += Math::conj(local(r, c));
        }
      }

#ifdef RODIN_USE_OPENMP
      /**
       * @brief Assembles @p out with OpenMP.
       *
       * Each thread owns a contiguous range of cells and scatters into the
       * shared value array; the additions that collide on a shared degree of
       * freedom are made atomically.
       *
       * @tparam KernelType Cell kernel, constructible from the two spaces.
       */
      template <class KernelType, class IterationType, class TrialFES, class TestFES>
      void assemble(MatrixType& out, const TrialFES& trialFES, const TestFES& testFES,
        const IterationType& seq, size_t d, Index count, int threadCount) const
      {
        if (!isValid(out, trialFES, testFES, count))
        {
          const size_t capacity = getCapacity(trialFES, testFES);
          std::vector<std::vector<Eigen::Triplet<ScalarType>>> chunks(
            static_cast<size_t>(threadCount));

#pragma omp parallel num_threads(threadCount)
          {
            KernelType kernel(trialFES, testFES);
            Math::Matrix<ScalarType> local;
            auto& triplets = chunks[static_cast<size_t>(omp_get_thread_num())];
            triplets.reserve(capacity / std::max(1, threadCount));
#pragma omp for
            for (Index i = 0; i < count; ++i)
            {
              const auto polytope = seq.getPolytope(i);
              kernel.compute(local, polytope);
              emplace(triplets, local, trialFES, testFES, d, i);
            }
          }

          std::vector<Eigen::Triplet<ScalarType>> triplets;
          triplets.reserve(capacity);
          for (auto& chunk : chunks)
            triplets.insert(triplets.end(), std::make_move_iterator(chunk.begin()),
              std::make_move_iterator(chunk.end()));
          setFromTriplets(out, triplets, trialFES, testFES);
          build(out, trialFES, testFES, d, count);
          return;
        }

        auto* const values = out.valuePtr();
        const Index nonZeroCount = out.nonZeros();
#pragma omp parallel num_threads(threadCount)
        {
#pragma omp for
          for (Index i = 0; i < nonZeroCount; ++i)
            values[i] = ScalarType(0);

          KernelType kernel(trialFES, testFES);
          Math::Matrix<ScalarType> local;
#pragma omp for
          for (Index i = 0; i < count; ++i)
          {
            const auto polytope = seq.getPolytope(i);
            kernel.compute(local, polytope);
            const auto& rows = testFES.getDOFs(d, i);
            const auto& cols = trialFES.getDOFs(d, i);
            size_t k = m_offsets[i];
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

      template <class TrialFES, class TestFES>
      static size_t getCapacity(const TrialFES& trialFES, const TestFES& testFES)
      {
        return testFES.getSize() * std::log(trialFES.getSize());
      }

      template <class TrialFES, class TestFES>
      static void emplace(std::vector<Eigen::Triplet<ScalarType>>& triplets,
        const Math::Matrix<ScalarType>& local, const TrialFES& trialFES,
        const TestFES& testFES, size_t d, Index i)
      {
        const auto& rows = testFES.getDOFs(d, i);
        const auto& cols = trialFES.getDOFs(d, i);
        for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
          for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
            triplets.emplace_back(rows(r), cols(c), Math::conj(local(r, c)));
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

      template <class TrialFES, class TestFES>
      bool isValid(const MatrixType& out, const TrialFES& trialFES,
        const TestFES& testFES, Index count) const
      {
        return m_built && out.isCompressed() &&
          static_cast<size_t>(out.rows()) == testFES.getSize() &&
          static_cast<size_t>(out.cols()) == trialFES.getSize() &&
          static_cast<size_t>(out.nonZeros()) == m_nonZeroCount &&
          m_offsets.size() == static_cast<size_t>(count) + 1;
      }

      /**
       * @brief Locates every cell-local entry in the value array of @p out.
       */
      template <class TrialFES, class TestFES>
      void build(const MatrixType& out, const TrialFES& trialFES, const TestFES& testFES,
        size_t d, Index count) const
      {
        const auto* const outer = out.outerIndexPtr();
        const auto* const inner = out.innerIndexPtr();
        m_offsets.resize(static_cast<size_t>(count) + 1);
        m_offsets[0] = 0;
        for (Index i = 0; i < count; ++i)
        {
          m_offsets[i + 1] = m_offsets[i] +
            static_cast<size_t>(testFES.getDOFs(d, i).size()) *
              static_cast<size_t>(trialFES.getDOFs(d, i).size());
        }
        m_indices.resize(m_offsets.back());
        for (Index i = 0; i < count; ++i)
        {
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
        }
        m_nonZeroCount = static_cast<size_t>(out.nonZeros());
        m_built = true;
      }

      mutable bool m_built = false;
      mutable size_t m_nonZeroCount = 0;
      mutable std::vector<size_t> m_offsets;
      mutable std::vector<Index> m_indices;
  };
}

#endif
