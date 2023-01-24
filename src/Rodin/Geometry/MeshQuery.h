/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MESHQUERY_H
#define RODIN_GEOMETRY_MESHQUERY_H

#include <functional>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  class MeshSelect
  {
    public:
      MeshSelect(const MeshBase& mesh)
        : m_mesh(mesh)
      {}

      MeshSelect& where(Attribute attr);

      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

    private:
      std::reference_wrapper<const MeshBase> m_mesh;
  };
}

#endif

