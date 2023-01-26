/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CONNECTIVITY_H
#define RODIN_GEOMETRY_CONNECTIVITY_H

#include <vector>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents the set of incidence relations.
   *
   * Stores the set of incidence relations:
   * @f[
   *  d \rightarrow d'
   * @f]
   * for a fixed pair of topological dimensions @f$ (d, d') @f$.
   */
  class Connectivity
  {
    public:
      Connectivity(size_t d, size_t dp)
        : m_left(d), m_right(dp), m_offsets{0}
      {}

      size_t getLeft() const
      {
        return m_left;
      }

      size_t getRight() const
      {
        return m_right;
      }

      Connectivity& connect(const std::vector<Index>& incidence)
      {
        m_offsets.push_back(incidence.size());
        for (const auto& idx : incidence)
          m_indices.push_back(idx);
        return *this;
      }

      Connectivity& setIndices(std::vector<Index> indices)
      {
        m_indices = std::move(indices);
        return *this;
      }

      Connectivity& setOffsets(std::vector<size_t> offsets)
      {
        m_offsets = std::move(offsets);
        return *this;
      }

      /**
       * @brief Gets the indices of the simplices of dimension @f$ d' @f$,
       * incident to the simplex @f$ (d, i) @f$.
       */
      std::vector<Index> getIncidence(Index idx) const
      {
        assert(idx + 1 < m_offsets.size());
        std::vector<Index> res(m_offsets[idx + 1] - m_offsets[idx]);
        size_t j = 0;
        for (size_t i = m_offsets[idx]; i < m_offsets[idx + 1]; i++)
          res[j++] = m_indices[i];
        return res;
      }

      void build(size_t d)
      {
        assert(false);
      }

      void transpose(size_t d, size_t dp)
      {
        assert(false);
      }

      void intersection(size_t d, size_t dp)
      {
        assert(false);
      }

    private:
      size_t m_left, m_right;
      std::vector<Index> m_indices;
      std::vector<size_t> m_offsets;
  };
}

#endif
