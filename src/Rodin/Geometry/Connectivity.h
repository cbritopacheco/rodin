/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CONNECTIVITY_H
#define RODIN_GEOMETRY_CONNECTIVITY_H

#include <vector>

#include "Rodin/Array.h"

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
      Connectivity(size_t d, size_t dp, size_t n = 0)
        : m_left(d), m_right(dp), m_connectivity(n)
      {}

      Connectivity(const Connectivity&) = default;

      Connectivity(Connectivity&&) = default;

      size_t getLeft() const
      {
        return m_left;
      }

      size_t getRight() const
      {
        return m_right;
      }

      size_t getSize() const
      {
        return m_connectivity.size();
      }

      Connectivity& setSize(size_t size)
      {
        m_connectivity.resize(size);
        return *this;
      }

      Connectivity& connect(Index idx, const Array<Index>& incidence)
      {
        if (idx + 1 > m_connectivity.size())
          m_connectivity.resize(idx + 1);
        assert(idx < m_connectivity.size());
        m_connectivity[idx] = incidence;
        return *this;
      }

      /**
       * @brief Gets the indices of the simplices of dimension @f$ d' @f$,
       * incident to the simplex @f$ (d, i) @f$.
       */
      const Array<Index>& getIncidence(Index idx) const
      {
        assert(idx < getSize());
        return m_connectivity[idx];
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
      std::vector<Array<Index>> m_connectivity;
  };
}

#endif
