/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEXINDEXED_H
#define RODIN_GEOMETRY_SIMPLEXINDEXED_H

#include <deque>
#include <vector>
#include <unordered_map>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  template <class T>
  class SimplexIndexed
  {
    public:
      using Map = std::unordered_map<Index, T>;

      SimplexIndexed() = default;

      SimplexIndexed(size_t meshDimension)
        : m_tracked(meshDimension + 1)
      {}

      SimplexIndexed(const SimplexIndexed&) = default;

      SimplexIndexed(SimplexIndexed&&) = default;

      SimplexIndexed& operator=(const SimplexIndexed&) = default;

      SimplexIndexed& operator=(SimplexIndexed&&) = default;

      auto begin(size_t d)
      {
        return m_tracked[d].begin();
      }

      auto end(size_t d)
      {
        return m_tracked[d].end();
      }

      auto begin(size_t d) const
      {
        return m_tracked[d].begin();
      }

      auto end(size_t d) const
      {
        return m_tracked[d].end();
      }

      auto cbegin(size_t d) const
      {
        return m_tracked[d].begin();
      }

      auto cend(size_t d) const
      {
        return m_tracked[d].end();
      }

      inline
      SimplexIndexed& initialize(size_t meshDimension)
      {
        m_tracked.resize(meshDimension + 1);
        return *this;
      }

      inline
      SimplexIndexed& track(size_t d, Index i, T&& value)
      {
        assert(m_tracked.size() > 0);
        assert(d < m_tracked.size());
        m_tracked[d][i] = std::move(value);
        return *this;
      }

      inline
      SimplexIndexed& track(size_t d, Index i, const T& value)
      {
        assert(m_tracked.size() > 0);
        assert(d < m_tracked.size());
        m_tracked[d][i] = value;
        return *this;
      }

      inline
      T& at(size_t d, Index i)
      {
        assert(m_tracked.size() > 0);
        assert(isTracked(d, i));
        return m_tracked[d].at(i);
      }

      inline
      const T& at(size_t d, Index i) const
      {
        assert(m_tracked.size() > 0);
        assert(d < m_tracked.size());
        assert(isTracked(d, i));
        return m_tracked[d].at(i);
      }

      inline
      bool isTracked(size_t d, Index i) const
      {
        assert(m_tracked.size() > 0);
        assert(d < m_tracked.size());
        return m_tracked[d].count(i) > 0;
      }

      inline
      void reserve(size_t d, size_t count)
      {
        assert(m_tracked.size() > 0);
        assert(d < m_tracked.size());
        m_tracked[d].reserve(count);
      }

      inline
      void clear()
      {
        m_tracked.clear();
      }

      inline
      void clear(size_t d)
      {
        m_tracked[d].clear();
      }

    private:
      std::vector<std::unordered_map<Index, T>> m_tracked;
  };
}

#endif

