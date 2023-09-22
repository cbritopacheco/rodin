/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GEOMETRYMAP_H
#define RODIN_VARIATIONAL_GEOMETRYMAP_H

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "ForwardDecls.h"
#include "FiniteElement.h"
#include "FiniteElementSpace.h"


namespace Rodin::Variational
{
  template <class T>
  class GeometryMap
  {
    using GeometryType = Geometry::Polytope::Geometry;
    public:
      GeometryMap(std::initializer_list<std::pair<GeometryType, T>> l)
      {
        for (const auto& v : l)
        {
          const size_t g = static_cast<size_t>(v.first);
          assert(g > 0);
          assert(g < m_map.size());
          m_map[g] = v.second;
        }
      }

      T& operator[](GeometryType geom)
      {
        const size_t g = static_cast<size_t>(geom);
        assert(g > 0);
        assert(g < m_map.size());
        return m_map[g];
      }

      const T& operator[](GeometryType geom) const
      {
        const size_t g = static_cast<size_t>(geom);
        assert(g > 0);
        assert(g < m_map.size());
        return m_map[g];
      }

    private:
      std::array<T, Geometry::Polytope::Geometries.size()> m_map;
  };
}


#endif
