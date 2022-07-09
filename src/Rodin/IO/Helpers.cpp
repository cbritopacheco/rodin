#include "Helpers.h"

#include <iostream>

namespace Rodin::IO
{
   std::map<std::string, IO::MeshFormat> getMeshFileHeaders()
   {
      std::map<std::string, IO::MeshFormat> res;
      res.insert({"MFEM mesh v1.0",          IO::MeshFormat::MFEM});
      res.insert({"MFEM mesh v1.2",          IO::MeshFormat::MFEM});
      res.insert({"MFEM NC mesh v1.0",       IO::MeshFormat::MFEM});
      res.insert({"MFEM mesh v1.1",          IO::MeshFormat::MFEM});
      res.insert({"$IO::MeshFormat",         IO::MeshFormat::GMSH});
      res.insert({"MeshVersionFormatted 1",  IO::MeshFormat::MEDIT});
      res.insert({"MeshVersionFormatted 2",  IO::MeshFormat::MEDIT});
      return res;
   }

   std::optional<IO::MeshFormat> getMeshFormat(std::istream& input, bool seekBeg)
   {
      assert(input);

      std::string meshType;
      input >> std::ws;
      std::getline(input, meshType);
      if (seekBeg)
      {
         input.clear();
         input.seekg(0, std::ios::beg);
      }

      // Check for, and remove, a trailing '\\r' from and std::string.
      if (!meshType.empty() && *meshType.rbegin() == '\r')
         meshType.resize(meshType.size() - 1);

      if (MeshFileHeaders.count(meshType))
      {
         return MeshFileHeaders.at(meshType);
      }
      else
      {
         return {};
      }
   }
}

namespace Rodin::IO::Medit
{
   boost::bimap<std::string, SolKeyword> getSolKeywordMap()
   {
      boost::bimap<std::string, SolKeyword> res;
      res.insert({"SolAtVertices",          SolKeyword::SolAtVertices});
      res.insert({"SolAtEdges",             SolKeyword::SolAtEdges});
      res.insert({"SolAtTriangles",         SolKeyword::SolAtTriangles});
      res.insert({"SolAtQuadrilaterals",    SolKeyword::SolAtQuadrilaterals});
      res.insert({"SolAtTetrahedra",        SolKeyword::SolAtTetrahedra});
      res.insert({"SolAtPentahedra",        SolKeyword::SolAtPentahedra});
      res.insert({"SolAtHexahedra",         SolKeyword::SolAtHexahedra});
      return res;
   }

   boost::bimap<std::string, EntityKeyword> getEntityKeywordMap()
   {
      boost::bimap<std::string, EntityKeyword> res;
      res.insert({"Vertices",          EntityKeyword::Vertices});
      res.insert({"Edges",             EntityKeyword::Edges});
      res.insert({"Triangles",         EntityKeyword::Triangles});
      res.insert({"SolAtTetrahedra",   EntityKeyword::Tetrahedra});
      return res;
   }
}
