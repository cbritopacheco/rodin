/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/algorithm/string.hpp>

#include "Rodin/Alert.h"

#include "MeshLoader.h"

namespace Rodin::IO
{
   void MeshLoader<MeshFormat::MFEM, Traits::Serial>::load(std::istream& is)
   {
      assert(is);
      auto fmt = getMeshFormat(is);
      is.clear();
      is.seekg(0, std::ios::beg);
      if (!fmt.has_value())
         Alert::Exception("Unrecognized mesh format.").raise();
      if (fmt == MeshFormat::MFEM)
      {
         getObject() = Rodin::Mesh<Traits::Serial>(mfem::Mesh(is, 0, 1, getFixOrientation()));
      }
      else
      {
         Alert::Exception("Cannot determine MFEM format version.").raise();
      }
   }

   void MeshLoader<MeshFormat::GMSH, Traits::Serial>::load(std::istream& is)
   {
      assert(is);
      auto fmt = getMeshFormat(is);
      is.clear();
      is.seekg(0, std::ios::beg);
      if (!fmt.has_value())
         Alert::Exception("Unrecognized mesh format.").raise();
      if (fmt == MeshFormat::GMSH)
      {
         getObject() = Rodin::Mesh<Traits::Serial>(mfem::Mesh(is, 0, 1, getFixOrientation()));
      }
      else
      {
         Alert::Exception("Cannot determine GMSH format version.").raise();
      }
   }

   void MeshLoader<MeshFormat::MEDIT, Traits::Serial>::load(std::istream& is)
   {
      std::string line;
      while (std::getline(is, line))
      {
         boost::algorithm::trim(line);
         if (!line.empty())
            break;
      }

      std::istringstream iss(line);
      std::string kw;
      iss >> kw;

      int version = 0;
      if (Medit::KeywordMap.left.count(kw))
      {
         switch (Medit::KeywordMap.left.at(kw))
         {
            case Medit::Keyword::MeshVersionFormatted:
            {
               if (!(iss >> version)) // Version is not on this line
               {
                  // Try next line
                  std::getline(is, line);
                  if (!(std::istringstream(line) >> version))
                  {
                     Alert::Exception()
                        << "Bad mesh format. Encountered \""
                        << line << "\" after keyword \""
                        << Medit::Keyword::MeshVersionFormatted
                        << "\"."
                        << Alert::Raise;
                  }
               }
               break;
            }
            default:
            {
               Alert::Exception()
                  << "Bad mesh format. Encountered keyword \""
                  << kw
                  << "\" before keyword \""
                  << Medit::Keyword::MeshVersionFormatted
                  << "\"."
                  << Alert::Raise;
               break;
            }
         }
      }
      else
      {
         Alert::Exception()
            << "Bad mesh format. Encountered \""
            << line << "\" before keyword \""
            << Medit::Keyword::MeshVersionFormatted
            << "\"."
            << Alert::Raise;
      }

      assert(version >= 2 && version <= 3);

      mfem::Mesh mfemMesh;

      std::map<Medit::Keyword, std::optional<std::istream::pos_type>> pos;
      std::map<Medit::Keyword, size_t> count;

      int spaceDim = 0;
      while (std::getline(is, line))
      {
         boost::algorithm::trim(line);
         if (line.empty())
            continue;

         std::istringstream iss(line);
         std::string kw;
         iss >> kw;
         if (kw == Medit::KeywordMap.right.at(Medit::Keyword::Dimension))
         {
            if (!(iss >> spaceDim)) // Dimension not on this line
            {
               // Try next line
               std::getline(is, line);
               if (!(std::istringstream(line) >> spaceDim))
               {
                  Alert::Exception()
                     << "Bad solution format. Encountered \""
                     << line << "\" after keyword \""
                     << Medit::Keyword::Dimension
                     << "\"."
                     << Alert::Raise;
               }
            }
            if (spaceDim < 2 || spaceDim > 3)
               Alert::Exception() << "Invalid dimension " << spaceDim << Alert::Raise;
            break;
         }
         else
         {
            Alert::Exception()
               << "Bad solution format. Encountered keyword \"" << kw
               << "\" before keyword \""
               << Medit::Keyword::Dimension
               << "\"."
               << Alert::Raise;
            return;
         }
      }

      while (std::getline(is, line))
      {
         boost::algorithm::trim(line);
         if (line.empty())
            continue;

         std::string kw;
         std::istringstream(line) >> kw;

         if (!std::isalpha(kw[0]))
            continue;

         if (Medit::KeywordMap.left.count(kw))
         {
            switch (Medit::KeywordMap.left.at(kw))
            {
               case Medit::Keyword::Tetrahedra:
               case Medit::Keyword::Triangles:
               case Medit::Keyword::Edges:
               case Medit::Keyword::Vertices:
               {
                  auto ent = Medit::KeywordMap.left.at(kw);
                  std::getline(is, line);
                  pos[ent] = is.tellg();
                  std::istringstream lss(line);
                  lss >> count[ent];
                  continue;
               }
               case Medit::Keyword::End:
               {
                  break;
               }
               default:
               {
                  Alert::Warning()
                     << "Ignoring unrecognized keyword \"" << kw << "\"."
                     << Alert::Raise;
                  continue;
               }
            }
            break;
         }
         else
         {
            Alert::Warning()
               << "Ignoring unrecognized keyword \"" << kw << "\"."
               << Alert::Raise;
            continue;
         }
      }

      // Set entities with zero count
      for (const auto& [_, kw] : Medit::KeywordMap)
      {
         switch (kw)
         {
            case Medit::Keyword::Tetrahedra:
            case Medit::Keyword::Triangles:
            case Medit::Keyword::Edges:
            case Medit::Keyword::Vertices:
            {
               if (!count.count(kw))
                  count[kw] = 0;
               break;
            }
            default:
            {
               // No-op
            }
         }
      }

      // Infer type of mesh
      bool isSurfaceMesh;
      if (spaceDim == 3 && count.at(Medit::Keyword::Tetrahedra) == 0)
      {
         // It's a surface mesh
         mfemMesh = mfem::Mesh(
               spaceDim - 1, // Dimension
               count.at(Medit::Keyword::Vertices),
               count.at(Medit::Keyword::Triangles),
               count.at(Medit::Keyword::Edges),
               spaceDim // Space dimension
               );
         isSurfaceMesh = true;
      }
      else if (spaceDim == 2 && count.at(Medit::Keyword::Triangles) == 0)
      {
         // It's a surface mesh
         mfemMesh = mfem::Mesh(
               spaceDim - 1, // Dimension
               count.at(Medit::Keyword::Vertices),
               count.at(Medit::Keyword::Edges),
               spaceDim // Space dimension
               );
         isSurfaceMesh = true;
      }
      else if (spaceDim == 3)
      {
         // It's a volume mesh
         mfemMesh = mfem::Mesh(
               spaceDim,
               count.at(Medit::Keyword::Vertices),
               count.at(Medit::Keyword::Tetrahedra),
               count.at(Medit::Keyword::Triangles),
               spaceDim
               );
         isSurfaceMesh = false;
      }
      else if (spaceDim == 2)
      {
         // It's a volume mesh
         mfemMesh = mfem::Mesh(
               spaceDim,
               count.at(Medit::Keyword::Vertices),
               count.at(Medit::Keyword::Triangles),
               count.at(Medit::Keyword::Edges),
               spaceDim
               );
         isSurfaceMesh = false;
      }
      else
      {
         Alert::Exception("Unhandled case.").raise();
      }

      is.clear();
      for (const auto& [_, kw] : Medit::KeywordMap)
      {
         if (!pos.count(kw))
            continue;
         auto g = pos.at(kw);
         if (!g.has_value())
            continue;
         is.seekg(*g);

         if (!count.count(kw))
         {
            Alert::Exception()
               << "Failed to parse \"" << kw
               << "\" count."
               << Alert::Raise;
         }

         switch (kw)
         {
            case Medit::Keyword::Vertices:
            {
               // Read all vertices
               for (size_t i = 0; i < count.at(kw); i++)
               {
                  if(!std::getline(is, line))
                     Alert::Exception("Bad mesh format.").raise();
                  std::istringstream lss(line);
                  double coords[3];
                  for (int j = 0; j < spaceDim; j++)
                     lss >> coords[j];
                  mfemMesh.AddVertex(coords[0], coords[1], coords[2]);
                  // We ignore the reference
               }
               break;
            }
            case Medit::Keyword::Triangles:
            {
               assert(spaceDim >= 2);
               // Read all triangles
               for (size_t i = 0; i < count.at(kw); i++)
               {
                  if(!std::getline(is, line))
                     Alert::Exception("Bad mesh format.").raise();
                  std::istringstream lss(line);
                  int v1, v2, v3, ref;
                  lss >> v1 >> v2 >> v3 >> ref;
                  switch (spaceDim)
                  {
                     case 2:
                     {
                        mfemMesh.AddTriangle(v1 - 1, v2 - 1, v3 - 1, ref);
                        break;
                     }
                     case 3:
                     {
                        if (isSurfaceMesh)
                           mfemMesh.AddTriangle(v1 - 1, v2 - 1, v3 - 1, ref);
                        else
                           mfemMesh.AddBdrTriangle(v1 - 1, v2 - 1, v3 - 1, ref);
                        break;
                     }
                     default:
                     {
                        Alert::Exception("Unhandled case.").raise();
                     }
                  }
               }
               break;
            }
            case Medit::Keyword::Tetrahedra:
            {
               assert(spaceDim >= 3);
               // Read all tetrahedra
               for (size_t i = 0; i < count.at(kw); i++)
               {
                  if(!std::getline(is, line))
                     Alert::Exception("Bad mesh format.").raise();
                  std::istringstream lss(line);
                  int v1, v2, v3, v4, ref;
                  lss >> v1 >> v2 >> v3 >> v4 >> ref;
                  mfemMesh.AddTet(v1 - 1, v2 - 1, v3 - 1, v4 - 1, ref);
               }
               break;
            }
            case Medit::Keyword::Edges:
            {
               assert(spaceDim >= 2);
               // Read all edges
               for (size_t i = 0; i < count.at(kw); i++)
               {
                  if(!std::getline(is, line))
                     Alert::Exception("Bad mesh format.").raise();
                  std::istringstream lss(line);
                  int v1, v2, ref;
                  lss >> v1 >> v2 >> ref;
                  switch (spaceDim)
                  {
                     case 2: // Planar mesh
                     case 3: // Surface mesh
                     {
                        if (isSurfaceMesh)
                           mfemMesh.AddBdrSegment(v1 - 1, v2 - 1, ref);
                        break;
                     }
                     default:
                     Alert::Exception("Unhandled case.").raise();
                  }
               }
               break;
            }
            default:
            {
               Alert::Warning()
                  << "Ignoring unrecognized entity \""
                  << kw
                  << " \"."
                  << Alert::Raise;
            }
         }
      }

      mfemMesh.FinalizeMesh(0, getFixOrientation());

      getObject() = Rodin::Mesh<Traits::Serial>(std::move(mfemMesh));
   }
}
