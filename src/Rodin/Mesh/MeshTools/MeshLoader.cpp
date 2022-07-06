/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "MeshLoader.h"

namespace Rodin::MeshTools
{
   const std::map<std::string, MeshFormat> MeshLoaderBase::FILE_HEADERS =
   {
      {"MFEM mesh v1.0",      MeshFormat::MFEM},
      {"MFEM mesh v1.2",      MeshFormat::MFEM},
      {"MFEM NC mesh v1.0",   MeshFormat::MFEM},
      {"MFEM mesh v1.1",      MeshFormat::MFEM},
      {"$MeshFormat",         MeshFormat::GMSH}
   };

   std::optional<MeshFormat> MeshLoaderBase::getMeshFormat(std::istream& input)
   {
      assert(input);

      std::string meshType;
      input >> std::ws;
      std::getline(input, meshType);
      input.clear();
      input.seekg(0, std::ios::beg);

      // Check for, and remove, a trailing '\\r' from and std::string.
      if (!meshType.empty() && *meshType.rbegin() == '\r')
         meshType.resize(meshType.size() - 1);

      if (MeshLoaderBase::FILE_HEADERS.count(meshType))
      {
         return MeshLoaderBase::FILE_HEADERS.at(meshType);
      }
      else
      {
         return {};
      }
   }

   IO::Status MeshLoader<MeshFormat::MFEM>::load(std::istream& is)
   {
      assert(is);

      SetEmpty();

      if (getMeshFormat(is) == MeshFormat::MFEM)
      {
         mfem::Mesh::Load(is, 0, 1, getFixOrientation());
         return {true, {}};
      }
      else
      {
         return {false, IO::Error{"Cannot determine MFEM mesh format version."}};
      }
   }

   IO::Status MeshLoader<MeshFormat::GMSH>::load(std::istream& is)
   {
      assert(is);

      SetEmpty();

      if (getMeshFormat(is) == MeshFormat::GMSH)
      {
         mfem::Mesh::Load(is, 0, 1, getFixOrientation());
         return {true, {}};
      }
      else
      {
         return {false, IO::Error{"Cannot determine Gmsh mesh format version."}};
      }
   }
}
