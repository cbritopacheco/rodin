/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHPRINTER_H
#define RODIN_EXTERNAL_MMG_MESHPRINTER_H

#include "Rodin/IO/Printer.h"

#include "Mesh.h"

namespace Rodin::External::MMG
{
  class MeshPrinter : public IO::Printer<MMG::Mesh>
  {
    public:
       MeshPrinter(const MMG::Mesh& mesh)
          : m_mesh(mesh)
       {}

       void print(std::ostream& os) override;

       const MMG::Mesh& getObject() const override
       {
          return m_mesh;
       }

    private:
       const MMG::Mesh& m_mesh;
  };
}

#endif

