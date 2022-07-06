/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Loader.h"

namespace Rodin::MeshTools
{
   const std::set<std::string> LoaderBase::MFEM_HEADERS =
   {
      "MFEM mesh v1.0",
      "MFEM mesh v1.2",
      "MFEM NC mesh v1.0",
      "MFEM mesh v1.1"
   };

   const std::set<std::string> LoaderBase::GMSH_HEADERS =
   {
      "$MeshFormat"
   };

   std::optional<MeshFormat> LoaderBase::getMeshFormat(std::istream& input)
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

      if (LoaderBase::MFEM_HEADERS.count(meshType))
      {
         return MeshFormat::MFEM;
      }
      else
      {
         return {};
      }
   }

   LoaderBase::Status Loader<MeshFormat::MFEM>::load()
   {
      assert(getInput());

      SetEmpty();

      if (getMeshFormat(getInput()) == MeshFormat::MFEM)
      {
         mfem::Mesh::Load(getInput(), 0, 1, getFixOrientation());
         return {true, {}};
      }
      else
      {
         return {false, Error{"Cannot determine MFEM mesh format version."}};
      }
   }

   LoaderBase::Status Loader<MeshFormat::GMSH>::load()
   {
      assert(getInput());

      SetEmpty();

      if (getMeshFormat(getInput()) == MeshFormat::GMSH)
      {
         mfem::Mesh::Load(getInput(), 0, 1, getFixOrientation());
         return {true, {}};
      }
      else
      {
         return {false, Error{"Cannot determine Gmsh mesh format version."}};
      }
   }
}
