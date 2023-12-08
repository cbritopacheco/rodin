/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHLOADER_H
#define RODIN_EXTERNAL_MMG_MESHLOADER_H

#include "Rodin/IO/MEDIT.h"

#include "Mesh.h"

namespace Rodin::External::MMG
{
  class MeshLoader : public IO::MeshLoader<IO::FileFormat::MEDIT, Context::Sequential>
  {
   public:
     using Parent = IO::MeshLoader<IO::FileFormat::MEDIT, Context::Sequential>;

     MeshLoader(MMG::Mesh& mesh)
       : Parent(mesh),
         m_mesh(mesh)
     {}

     void load(std::istream& is) override;

     MMG::Mesh& getObject() override
     {
       return m_mesh.get();
     }

   private:
     std::reference_wrapper<MMG::Mesh> m_mesh;
  };
}

#endif
