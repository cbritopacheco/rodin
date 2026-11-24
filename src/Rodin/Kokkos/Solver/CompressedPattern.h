/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_COMPRESSEDPATTERN_H
#define RODIN_KOKKOS_SOLVER_COMPRESSEDPATTERN_H

/**
 * @file
 * @brief Symbolic compressed patterns and their factor constructions.
 */

#include <algorithm>
#include <vector>
#include "Rodin/Math/SparseMatrix.h"

namespace Rodin::Kokkos::Solver
{
  /**
   * @brief Compressed pattern of a square matrix, without values.
   *
   * ## Architecture
   *
   * Transpose retains the input positions of the transposed entries. Cholesky
   * constructs an elimination tree and its lower-factor pattern. LU propagates
   * symbolic row reach through previously eliminated columns. The constructors
   * operate on valid compressed arrays and permutations on the host; no numeric
   * values or device allocation participate.
   */
  struct CompressedPattern
  {
    /// @brief Start of every compressed line, followed by the entry count.
    std::vector<::Rodin::Math::SparseIndex> starts;

    /// @brief Index of every entry in its compressed line.
    std::vector<::Rodin::Math::SparseIndex> indices;

    /**
     * @brief Transposes the pattern of a square matrix.
     *
     * The indices of every line of the result are sorted.
     *
     * @param n Order of the matrix.
     * @param starts Line starts of the compressed arrays.
     * @param indices Indices of the compressed arrays.
     * @param[out] positions Position in @p indices of every entry of the result.
     * @returns Pattern of the transpose.
     */
    static CompressedPattern Transpose(Index n,
      const ::Rodin::Math::SparseIndex* starts,
      const ::Rodin::Math::SparseIndex* indices,
      std::vector<::Rodin::Math::SparseIndex>& positions)
    {
      const Index nonzeros = static_cast<Index>(starts[n]);
      CompressedPattern res;
      res.starts.assign(n + 1, 0);
      res.indices.resize(nonzeros);
      positions.resize(nonzeros);
      for (Index k = 0; k < nonzeros; ++k)
        ++res.starts[indices[k] + 1];
      for (Index i = 0; i < n; ++i)
        res.starts[i + 1] += res.starts[i];
      std::vector<::Rodin::Math::SparseIndex> next(
        res.starts.begin(), res.starts.end() - 1);
      for (Index line = 0; line < n; ++line)
      {
        for (::Rodin::Math::SparseIndex k = starts[line]; k < starts[line + 1]; ++k)
        {
          const ::Rodin::Math::SparseIndex at = next[indices[k]]++;
          res.indices[at] = static_cast<::Rodin::Math::SparseIndex>(line);
          positions[at] = k;
        }
      }
      return res;
    }

    /**
     * @brief Computes the CSR pattern of the Cholesky factor of
     * @f$ Q^T (A + A^T) Q @f$.
     *
     * The elimination tree is built first, and the pattern of every row is
     * then the set of nodes reached by climbing the tree from the entries of
     * the row below the diagonal. The cost is proportional to the number of
     * entries of the factor. The diagonal closes every row.
     *
     * @param n Order of the matrix.
     * @param starts Line starts of the compressed arrays of @f$ A @f$ or @f$ A^T @f$.
     * @param indices Indices of the compressed arrays.
     * @param q Original index of every position of the ordering.
     * @returns Pattern of the lower triangular factor, rows sorted.
     */
    static CompressedPattern Cholesky(Index n, const ::Rodin::Math::SparseIndex* starts,
      const ::Rodin::Math::SparseIndex* indices,
      const std::vector<::Rodin::Math::SparseIndex>& q)
    {
      std::vector<::Rodin::Math::SparseIndex> position(n);
      for (Index k = 0; k < n; ++k)
        position[q[k]] = static_cast<::Rodin::Math::SparseIndex>(k);

      // Strictly lower rows of the permuted symmetric pattern.
      CompressedPattern lower;
      lower.starts.assign(n + 1, 0);
      for (Index line = 0; line < n; ++line)
      {
        for (::Rodin::Math::SparseIndex k = starts[line]; k < starts[line + 1]; ++k)
        {
          const ::Rodin::Math::SparseIndex a = position[line];
          const ::Rodin::Math::SparseIndex b = position[indices[k]];
          if (a != b)
            ++lower.starts[std::max(a, b) + 1];
        }
      }
      for (Index i = 0; i < n; ++i)
        lower.starts[i + 1] += lower.starts[i];
      lower.indices.resize(lower.starts[n]);
      std::vector<::Rodin::Math::SparseIndex> next(
        lower.starts.begin(), lower.starts.end() - 1);
      for (Index line = 0; line < n; ++line)
      {
        for (::Rodin::Math::SparseIndex k = starts[line]; k < starts[line + 1]; ++k)
        {
          const ::Rodin::Math::SparseIndex a = position[line];
          const ::Rodin::Math::SparseIndex b = position[indices[k]];
          if (a != b)
            lower.indices[next[std::max(a, b)]++] = std::min(a, b);
        }
      }

      // Elimination tree, with path compression through the ancestors.
      std::vector<::Rodin::Math::SparseIndex> parent(n, -1);
      std::vector<::Rodin::Math::SparseIndex> ancestor(n, -1);
      for (Index row = 0; row < n; ++row)
      {
        for (::Rodin::Math::SparseIndex k = lower.starts[row];
             k < lower.starts[row + 1]; ++k)
        {
          ::Rodin::Math::SparseIndex node = lower.indices[k];
          while (node != -1 && node < static_cast<::Rodin::Math::SparseIndex>(row))
          {
            const ::Rodin::Math::SparseIndex up = ancestor[node];
            ancestor[node] = static_cast<::Rodin::Math::SparseIndex>(row);
            if (up == -1)
              parent[node] = static_cast<::Rodin::Math::SparseIndex>(row);
            node = up;
          }
        }
      }

      // Row patterns, as the subtrees reached from the entries of each row.
      CompressedPattern res;
      res.starts.reserve(n + 1);
      res.starts.push_back(0);
      std::vector<::Rodin::Math::SparseIndex> mark(n, -1);
      for (Index row = 0; row < n; ++row)
      {
        const ::Rodin::Math::SparseIndex r =
          static_cast<::Rodin::Math::SparseIndex>(row);
        const size_t first = res.indices.size();
        mark[row] = r;
        for (::Rodin::Math::SparseIndex k = lower.starts[row];
             k < lower.starts[row + 1]; ++k)
        {
          for (::Rodin::Math::SparseIndex node = lower.indices[k];
               node != -1 && mark[node] != r; node = parent[node])
          {
            mark[node] = r;
            res.indices.push_back(node);
          }
        }
        std::sort(res.indices.begin() + first, res.indices.end());
        res.indices.push_back(r);
        res.starts.push_back(
          static_cast<::Rodin::Math::SparseIndex>(res.indices.size()));
      }
      return res;
    }

    /**
     * @brief Computes the CSR pattern of @f$ (L - I) + U @f$ for the LU
     * factorization of @f$ PAQ @f$ without further pivoting.
     *
     * Every column of the factors is the set of rows reached in the graph of
     * the previous columns of @f$ L @f$ from the entries of the column of
     * @f$ PAQ @f$, as in the left-looking algorithm of Gilbert and Peierls.
     * The cost is that of the reachability searches, bounded by the number of
     * operations of the factorization. The diagonal is always kept.
     *
     * @param n Order of the matrix.
     * @param starts Column starts of the CSC arrays of @f$ A @f$.
     * @param indices Row indices of the CSC arrays of @f$ A @f$.
     * @param p Original row of every pivot row.
     * @param q Original column of every pivot column.
     * @returns Pattern of the bundled factors, rows sorted.
     */
    static CompressedPattern LU(Index n, const ::Rodin::Math::SparseIndex* starts,
      const ::Rodin::Math::SparseIndex* indices,
      const std::vector<::Rodin::Math::SparseIndex>& p,
      const std::vector<::Rodin::Math::SparseIndex>& q)
    {
      std::vector<::Rodin::Math::SparseIndex> pivot(n);
      for (Index k = 0; k < n; ++k)
        pivot[p[k]] = static_cast<::Rodin::Math::SparseIndex>(k);

      // Strictly lower part of the columns of L computed so far.
      CompressedPattern l;
      l.starts.reserve(n + 1);
      l.starts.push_back(0);

      // Every column of the bundle, rows unsorted.
      CompressedPattern t;
      t.starts.reserve(n + 1);
      t.starts.push_back(0);

      std::vector<::Rodin::Math::SparseIndex> mark(n, -1);
      std::vector<::Rodin::Math::SparseIndex> stack;
      for (Index column = 0; column < n; ++column)
      {
        const ::Rodin::Math::SparseIndex j =
          static_cast<::Rodin::Math::SparseIndex>(column);
        mark[j] = j;
        t.indices.push_back(j);
        for (::Rodin::Math::SparseIndex k = starts[q[column]];
             k < starts[q[column] + 1]; ++k)
        {
          const ::Rodin::Math::SparseIndex row = pivot[indices[k]];
          if (mark[row] != j)
          {
            mark[row] = j;
            stack.push_back(row);
          }
        }
        while (!stack.empty())
        {
          const ::Rodin::Math::SparseIndex node = stack.back();
          stack.pop_back();
          t.indices.push_back(node);
          if (node < j)
          {
            for (::Rodin::Math::SparseIndex k = l.starts[node]; k < l.starts[node + 1];
                 ++k)
            {
              const ::Rodin::Math::SparseIndex row = l.indices[k];
              if (mark[row] != j)
              {
                mark[row] = j;
                stack.push_back(row);
              }
            }
          }
        }
        for (size_t k = t.starts.back(); k < t.indices.size(); ++k)
        {
          if (t.indices[k] > j)
            l.indices.push_back(t.indices[k]);
        }
        l.starts.push_back(static_cast<::Rodin::Math::SparseIndex>(l.indices.size()));
        t.starts.push_back(static_cast<::Rodin::Math::SparseIndex>(t.indices.size()));
      }

      // The columns are visited in order, so the transpose has sorted rows.
      std::vector<::Rodin::Math::SparseIndex> positions;
      return Transpose(n, t.starts.data(), t.indices.data(), positions);
    }
  };
}

#endif
