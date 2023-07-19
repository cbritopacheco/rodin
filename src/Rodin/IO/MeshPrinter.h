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
  template <class Trait>
  class MeshPrinterBase : public IO::Printer<Rodin::Geometry::Mesh<Trait>>
  {
    public:
      MeshPrinterBase(const Geometry::Mesh<Context::Serial>& mesh)
        : m_mesh(mesh)
      {}

      const Rodin::Geometry::Mesh<Context::Serial>& getObject() const override
      {
        return m_mesh;
      }

    private:
      const Rodin::Geometry::Mesh<Context::Serial>& m_mesh;
  };
}

#endif

