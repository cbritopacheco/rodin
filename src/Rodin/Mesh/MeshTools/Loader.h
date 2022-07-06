/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHTOOLS_LOADER_H
#define RODIN_MESH_MESHTOOLS_LOADER_H

#include <optional>
#include <mfem.hpp>
#include <boost/filesystem.hpp>

#include "ForwardDecls.h"

namespace Rodin::MeshTools
{
   class LoaderBase : protected mfem::Mesh
   {
      public:
         static const std::set<std::string> MFEM_HEADERS;

         static const std::set<std::string> GMSH_HEADERS;

         static std::optional<MeshFormat> getMeshFormat(std::istream& input);

         struct Error
         {
            std::string message;
         };

         struct Status
         {
            bool success;
            std::optional<Error> error;
         };

         LoaderBase(std::istream& is)
            : m_is(is),
              m_fixOrientation(true)
         {}

         std::istream& getInput() const
         {
            return m_is;
         }

         mfem::Mesh& getMesh()
         {
            return *this;
         }

         LoaderBase& fixOrientation(bool fixOrientation)
         {
            m_fixOrientation = fixOrientation;
            return *this;
         }

         virtual Status load() = 0;

      protected:
         bool getFixOrientation() const
         {
            return m_fixOrientation;
         }

      private:
         std::istream& m_is;
         bool m_fixOrientation;
   };

   using LoaderStatus = LoaderBase::Status;

   template <>
   class Loader<MeshFormat::MFEM> : public LoaderBase
   {
      public:
         Loader(std::istream& is)
            : LoaderBase(is)
         {}

         Loader(const Loader& other)
            :  LoaderBase(other.getInput()),
               m_parseTag(other.m_parseTag)
         {}

         Loader(Loader&& other)
            : LoaderBase(std::move(other)),
              m_parseTag(std::move(other.m_parseTag))
         {}

         Status load() override;

      private:
         std::string m_parseTag;
   };

   template <>
   class Loader<MeshFormat::GMSH> : public LoaderBase
   {
      public:
         Loader(std::istream& is)
            : LoaderBase(is)
         {}

         Loader(const Loader& other)
            :  LoaderBase(other.getInput())
         {}

         Loader(Loader&& other)
            : LoaderBase(std::move(other))
         {}

         Status load() override;
   };
}

#endif
