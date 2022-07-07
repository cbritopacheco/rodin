/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/algorithm/string.hpp>

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
      {"$MeshFormat",         MeshFormat::GMSH},
      {"MeshVersionFormatted 1", MeshFormat::MEDIT},
      {"MeshVersionFormatted 2", MeshFormat::MEDIT}
   };

   const std::map<std::string, MeshLoader<MeshFormat::MEDIT>::EntityKeyword>
      MeshLoader<MeshFormat::MEDIT>::s_entities =
   {
      {"Vertices",      EntityKeyword::Vertices},
      {"Triangles",     EntityKeyword::Triangles},
      {"Tetrahedra",    EntityKeyword::Tetrahedra},
      {"Edges",         EntityKeyword::Edges}
   };

   std::optional<MeshFormat> MeshLoaderBase::getMeshFormat(std::istream& input, bool seekBeg)
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

   IO::Status MeshLoader<MeshFormat::MEDIT>::load(std::istream& is)
   {
      SetEmpty();

      int spaceDim = 0;
      if (getMeshFormat(is, false) == MeshFormat::MEDIT)
      {
         std::string line;

         std::map<EntityKeyword, std::optional<std::istream::pos_type>> pos;
         std::map<EntityKeyword, size_t> count;

         while (std::getline(is, line))
         {
            boost::algorithm::trim(line);
            if (line.empty())
               continue;

            std::istringstream iss(line);
            std::string kw;
            iss >> kw;
            if (kw == "Dimension")
            {
               iss >> spaceDim;
               if (spaceDim < 2 || spaceDim > 3)
                  return {false, IO::Error{"Invalid mesh dimension: " + std::to_string(spaceDim)}};
            }
            else if (s_entities.count(kw))
            {
               auto ent = s_entities.at(kw);
               std::getline(is, line);
               pos[ent] = is.tellg();
               std::istringstream lss(line);
               lss >> count[ent];
            }
         }

         // Infer type of mesh
         bool isSurfaceMesh;
         if (spaceDim == 3 && (
                  count.count(EntityKeyword::Tetrahedra) == 0
                  || count.at(EntityKeyword::Tetrahedra) == 0))
         {
            // It's a surface mesh
            InitMesh(
                  spaceDim - 1, // Dimension
                  spaceDim, // Space dimension
                  count.at(EntityKeyword::Vertices),
                  count.at(EntityKeyword::Triangles),
                  count.at(EntityKeyword::Edges));
            isSurfaceMesh = true;
         }
         else if (spaceDim == 2 && (
                  count.count(EntityKeyword::Triangles) == 0
                  || count.at(EntityKeyword::Triangles) == 0))
         {
            // It's a surface mesh
            InitMesh(
                  spaceDim - 1, // Dimension
                  spaceDim, // Space dimension
                  count.at(EntityKeyword::Vertices),
                  count.at(EntityKeyword::Edges),
                  0);
            isSurfaceMesh = true;
         }
         else if (spaceDim == 3)
         {
            // It's a volume mesh
            InitMesh(
                  spaceDim,
                  spaceDim,
                  count.at(EntityKeyword::Vertices),
                  count.at(EntityKeyword::Tetrahedra),
                  count.at(EntityKeyword::Triangles)
                  );
            isSurfaceMesh = false;
         }
         else if (spaceDim == 2)
         {
            // It's a volume mesh
            InitMesh(
                  spaceDim,
                  spaceDim,
                  count.at(EntityKeyword::Vertices),
                  count.at(EntityKeyword::Triangles),
                  count.at(EntityKeyword::Edges)
                  );
            isSurfaceMesh = false;
         }
         else
         {
            Alert::Exception("Unhandled case.").raise();
         }

         is.clear();
         for (const auto& [kw, ent] : s_entities)
         {
            if (pos.count(ent))
            {
               auto g = pos.at(ent);
               if (g.has_value())
               {
                  is.seekg(*g);
                  if (count.count(ent))
                  {
                     switch (ent)
                     {
                        case EntityKeyword::Vertices:
                        {
                           // Read all vertices
                           for (size_t i = 0; i < count.at(ent); i++)
                           {
                              if(!std::getline(is, line))
                                 return IO::Status{false, IO::Error{"Bad mesh format."}};
                              std::istringstream lss(line);
                              double coords[3];
                              for (int j = 0; j < spaceDim; j++)
                                 lss >> coords[j];
                              AddVertex(coords[0], coords[1], coords[2]);
                              // We ignore the reference
                           }
                           break;
                        }
                        case EntityKeyword::Triangles:
                        {
                           assert(spaceDim >= 2);
                           // Read all triangles
                           for (size_t i = 0; i < count.at(ent); i++)
                           {
                              if(!std::getline(is, line))
                                 return IO::Status{false, IO::Error{"Bad mesh format."}};
                              std::istringstream lss(line);
                              int v1, v2, v3, ref;
                              lss >> v1 >> v2 >> v3 >> ref;
                              switch (spaceDim)
                              {
                                 case 2:
                                 {
                                    AddTriangle(v1 - 1, v2 - 1, v3 - 1, ref);
                                    break;
                                 }
                                 case 3:
                                 {
                                    if (isSurfaceMesh)
                                       AddTriangle(v1 - 1, v2 - 1, v3 - 1, ref);
                                    else
                                       AddBdrTriangle(v1 - 1, v2 - 1, v3 - 1, ref);
                                    break;
                                 }
                                 default:
                                 Alert::Exception("Unhandled case.").raise();
                              }
                           }
                           break;
                        }
                        case EntityKeyword::Tetrahedra:
                        {
                           assert(spaceDim >= 3);
                           // Read all tetrahedra
                           for (size_t i = 0; i < count.at(ent); i++)
                           {
                              if(!std::getline(is, line))
                                 return IO::Status{false, IO::Error{"Bad mesh format."}};
                              std::istringstream lss(line);
                              int v1, v2, v3, v4, ref;
                              lss >> v1 >> v2 >> v3 >> v4 >> ref;
                              AddTet(v1 - 1, v2 - 1, v3 - 1, v4 - 1, ref);
                           }
                           break;
                        }
                        case EntityKeyword::Edges:
                        {
                           assert(spaceDim >= 2);
                           // Read all edges
                           for (size_t i = 0; i < count.at(ent); i++)
                           {
                              if(!std::getline(is, line))
                                 return IO::Status{false, IO::Error{"Bad mesh format."}};
                              std::istringstream lss(line);
                              int v1, v2, ref;
                              lss >> v1 >> v2 >> ref;
                              switch (spaceDim)
                              {
                                 case 2: // Planar mesh
                                 case 3: // Surface mesh
                                 {
                                    AddBdrSegment(v1 - 1, v2 - 1, ref);
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
                           Alert::Exception("Unhandled Medit entity.").raise();
                        }
                     }
                  }
                  else
                  {
                     return {false, IO::Error{"Bad mesh formatting."}};
                  }
               }
            }
         }

         FinalizeMesh(0, getFixOrientation());

         return {true, {}};
      }
      else
      {
         return {false, IO::Error{"Cannot determine Medit mesh format version."}};
      }
   }
}
