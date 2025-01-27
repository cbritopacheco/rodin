/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_TRIANGULATION_H
#define RODIN_GEOMETRY_TRIANGULATION_H

#include "Polytope.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
  template <class T>
  class Simplexification;

  template <>
  class Simplexification<Polytope::Type>
  {
    public:
      constexpr
      Simplexification(Polytope::Type g)
        : m_geometry(g)
      {}

      constexpr
      Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      constexpr
      size_t getSize() const
      {
        return getSimplices().size();
      }

      virtual const std::vector<IndexArray>& getSimplices() const = 0;

    private:
      Polytope::Type m_geometry;
  };

  template <class MeshContext>
  class Simplexification<Mesh<MeshContext>>
  {
    public:
      using ContextType = MeshContext;
      using MeshType = Mesh<ContextType>;

      Simplexification(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      const MeshType& getMesh() const
      {
        return m_mesh.get();
      }

      virtual MeshType simplexify() = 0;

    private:
      std::reference_wrapper<const MeshType> m_mesh;
  };
}

#endif
