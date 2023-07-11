/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GEOMETRYINDEXED_H
#define RODIN_VARIATIONAL_GEOMETRYINDEXED_H

#include "Simplex.h"

namespace Rodin::Geometry
{
  template <class T>
  class GeometryIndexed
  {
    public:
      constexpr
      GeometryIndexed()
      {
        for (auto& v : m_map)
          v = T();
      }

      constexpr
      GeometryIndexed(std::initializer_list<std::pair<Polytope::Geometry, T>> l)
      {
        // All geometries must be handled exactly once
        assert(l.size() == Polytope::Geometries.size());

        for (const auto& v : l)
          m_map[static_cast<size_t>(v.first)] = v.second;
      }

      inline
      constexpr
      T& operator[](Polytope::Geometry geom)
      {
        const size_t g = static_cast<size_t>(geom);
        assert(g < m_map.size());
        return m_map[g];
      }

      inline
      constexpr
      const T& operator[](Polytope::Geometry geom) const
      {
        const size_t g = static_cast<size_t>(geom);
        assert(g < m_map.size());
        return m_map[g];
      }

      inline
      constexpr
      size_t size() const
      {
        return Polytope::Geometries.size();
      }

    private:
      std::array<T, Polytope::Geometries.size()> m_map;
  };
}


#endif

