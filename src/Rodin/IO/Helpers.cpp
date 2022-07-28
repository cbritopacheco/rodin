#include <boost/algorithm/string.hpp>

#include "Helpers.h"

namespace Rodin::IO
{
   std::ostream& operator<<(std::ostream& os, FileFormat fmt)
   {
      switch (fmt)
      {
         case FileFormat::MFEM:
         {
            os << "MFEM";
            break;
         }
         case FileFormat::MEDIT:
         {
            os << "MEDIT";
            break;
         }
         case FileFormat::GMSH:
         {
            os << "GMSH";
            break;
         }
      }
      return os;
   }

   std::map<std::string, IO::FileFormat> getMeshFileHeaders()
   {
      std::map<std::string, IO::FileFormat> res;
      res.insert({"MFEM mesh v1.0",          IO::FileFormat::MFEM});
      res.insert({"MFEM mesh v1.2",          IO::FileFormat::MFEM});
      res.insert({"MFEM NC mesh v1.0",       IO::FileFormat::MFEM});
      res.insert({"MFEM mesh v1.1",          IO::FileFormat::MFEM});
      res.insert({"$MeshFormat",             IO::FileFormat::GMSH});
      res.insert({"MeshVersionFormatted 1",  IO::FileFormat::MEDIT});
      res.insert({"MeshVersionFormatted 2",  IO::FileFormat::MEDIT});
      return res;
   }

   std::optional<IO::FileFormat> getMeshFormat(std::istream& input)
   {
      assert(input);

      std::string meshType;
      while (std::getline(input, meshType))
      {
         boost::algorithm::trim(meshType);
         if (!meshType.empty())
            break;
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
   boost::bimap<std::string, Keyword> getKeywordMap()
   {
      boost::bimap<std::string, Keyword> res;
      res.insert({"MeshVersionFormatted", Keyword::MeshVersionFormatted});
      res.insert({"Dimension",            Keyword::Dimension});
      res.insert({"Vertices",             Keyword::Vertices});
      res.insert({"Edges",                Keyword::Edges});
      res.insert({"Triangles",            Keyword::Triangles});
      res.insert({"Tetrahedra",           Keyword::Tetrahedra});
      res.insert({"SolAtVertices",        Keyword::SolAtVertices});
      res.insert({"SolAtEdges",           Keyword::SolAtEdges});
      res.insert({"SolAtTriangles",       Keyword::SolAtTriangles});
      res.insert({"SolAtQuadrilaterals",  Keyword::SolAtQuadrilaterals});
      res.insert({"SolAtTetrahedra",      Keyword::SolAtTetrahedra});
      res.insert({"SolAtPentahedra",      Keyword::SolAtPentahedra});
      res.insert({"SolAtHexahedra",       Keyword::SolAtHexahedra});
      res.insert({"End",                  Keyword::End});
      return res;
   }

   std::ostream& operator<<(std::ostream& os, Keyword kw)
   {
      assert(KeywordMap.right.count(kw));
      os << KeywordMap.right.at(kw);
      return os;
   }
}
