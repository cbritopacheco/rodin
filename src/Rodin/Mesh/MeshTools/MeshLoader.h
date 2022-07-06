/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHTOOLS_LOADER_H
#define RODIN_MESH_MESHTOOLS_LOADER_H

#include <map>
#include <optional>
#include <mfem.hpp>
#include <boost/filesystem.hpp>

#include "Rodin/IO/Loader.h"

#include "ForwardDecls.h"

namespace Rodin::MeshTools
{
   class MeshLoaderBase : public IO::Loader<mfem::Mesh>, protected mfem::Mesh
   {
      static const std::map<std::string, MeshFormat> FILE_HEADERS;

      public:
         static std::optional<MeshFormat> getMeshFormat(std::istream& input);

         MeshLoaderBase()
            : m_fixOrientation(true)
         {}

         mfem::Mesh& getMesh()
         {
            return *this;
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
         bool m_fixOrientation;
   };

   template <>
   class MeshLoader<MeshFormat::MFEM> : public MeshLoaderBase
   {
      public:
         IO::Status load(std::istream& is) override;

      private:
         std::string m_parseTag;
   };

   template <>
   class MeshLoader<MeshFormat::GMSH> : public MeshLoaderBase
   {
      public:
         IO::Status load(std::istream& is) override;
   };
}

#endif
