/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHLOADER_H
#define RODIN_MESH_MESHLOADER_H

#include <map>
#include <optional>
#include <mfem.hpp>
#include <boost/filesystem.hpp>

#include "Rodin/Mesh.h"
#include "Rodin/IO/Loader.h"

#include "ForwardDecls.h"
#include "Helpers.h"

namespace Rodin::IO
{
   template <class Trait>
   class MeshLoaderBase : public IO::Loader<Rodin::Mesh<Trait>>
   {
      public:
         MeshLoaderBase()
            : m_fixOrientation(true)
         {}

         Rodin::Mesh<Trait>& getObject() override
         {
            return m_mesh;
         }

         MeshLoaderBase& fixOrientation(bool fixOrientation)
         {
            m_fixOrientation = fixOrientation;
            return *this;
         }

      protected:
         bool getFixOrientation() const
         {
            return m_fixOrientation;
         }

      private:
         Rodin::Mesh<Traits::Serial> m_mesh;
         bool m_fixOrientation;
   };

   template <>
   class MeshLoader<IO::MeshFormat::MFEM, Traits::Serial>
      : public MeshLoaderBase<Traits::Serial>
   {
      public:
         IO::Status load(std::istream& is) override;

      private:
         std::string m_parseTag;
   };

   template <>
   class MeshLoader<IO::MeshFormat::GMSH, Traits::Serial>
      : public MeshLoaderBase<Traits::Serial>
   {
      public:
         IO::Status load(std::istream& is) override;
   };

   template <>
   class MeshLoader<IO::MeshFormat::MEDIT, Traits::Serial>
      : public MeshLoaderBase<Traits::Serial>
   {
      public:
         IO::Status load(std::istream& is) override;
   };
}

#endif
