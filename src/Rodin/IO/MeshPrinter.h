/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHPRINTER_H
#define RODIN_MESH_MESHPRINTER_H

#include <utility>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/IO/Printer.h"

#include "ForwardDecls.h"

namespace Rodin::IO
{
  template <class Context>
  class MeshPrinterBase : public IO::Printer<Geometry::Mesh<Context>>
  {
    public:
      using ContextType = Context;

      using ObjectType = Geometry::Mesh<ContextType>;

      MeshPrinterBase(const ObjectType& mesh)
        : m_mesh(mesh)
      {}

      const ObjectType& getObject() const override
      {
        return m_mesh.get();
      }

    private:
      std::reference_wrapper<const ObjectType> m_mesh;
  };
}

#endif

