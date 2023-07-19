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

#include "Types.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  template <class T>
  class PolytopeIndexed
  {
    public:
      PolytopeIndexed() = default;

      PolytopeIndexed(size_t meshDimension)
        : m_tracked(meshDimension + 1)
      {}

      PolytopeIndexed(const PolytopeIndexed&) = default;

      PolytopeIndexed(PolytopeIndexed&&) = default;

      PolytopeIndexed& operator=(const PolytopeIndexed&) = default;

      PolytopeIndexed& operator=(PolytopeIndexed&&) = default;

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

      size_t size(size_t d) const
      {
        return m_tracked[d].size();
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
      PolytopeIndexed& initialize(size_t meshDimension)
      {
        m_tracked.resize(meshDimension + 1);
        return *this;
      }

      inline
      PolytopeIndexed& track(const std::pair<size_t, Index>& p, T&& value)
      {
        const auto [d, i] = p;
        assert(m_tracked.size() > 0);
        assert(d < m_tracked.size());
        m_tracked[d][i] = std::move(value);
        return *this;
      }

      inline
      PolytopeIndexed& track(const std::pair<size_t, Index>& p, const T& value)
      {
        const auto [d, i] = p;
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
      auto find(size_t d, Index i) const
      {
        return m_tracked[d].find(i);
      }

      template <class Iterator>
      inline
      auto insert(Iterator&& it, const std::pair<size_t, Index>& p, const T& value)
      {
        auto [d, i] = p;
        return m_tracked[d].insert(std::forward<Iterator>(it), { i, value });
      }

      template <class Iterator>
      inline
      auto insert(Iterator&& it, const std::pair<size_t, Index>& p, T&& value)
      {
        auto [d, i] = p;
        return m_tracked[d].insert(std::forward<Iterator>(it), { i, std::move(value) });
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
        for (auto& m : m_tracked)
          m.clear();
      }

      inline
      void clear(size_t d)
      {
        m_tracked[d].clear();
      }

    private:
      std::vector<UnorderedMap<Index, T>> m_tracked;
  };
}

#endif

