/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHLOADER_H
#define RODIN_MESH_MESHLOADER_H

#include <utility>
#include <fstream>

#include "Rodin/IO/Loader.h"
#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @defgroup MeshLoaderSpecializations MeshLoader Template Specializations
   * @brief Template specializations of the MeshLoader class.
   * @see MeshLoader
   */

  /**
   * @brief Base class for mesh loader objects.
   */
  template <class Context>
  class MeshLoaderBase : public IO::Loader<Rodin::Geometry::Mesh<Context>>
  {
    public:
      using ContextType = Context;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = IO::Loader<ObjectType>;

      MeshLoaderBase(ObjectType& mesh)
        : m_mesh(mesh)
      {}

    protected:
      ObjectType& getObject() override
      {
        return m_mesh.get();
      }

    private:
      std::reference_wrapper<ObjectType> m_mesh;
  };
}

#endif
