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
   IO::Status MeshLoader<MeshFormat::MFEM, Traits::Serial>::load(std::istream& is)
   {
      assert(is);
      auto fmt = getMeshFormat(is);
      if (!fmt.has_value())
         Alert::Exception("Unrecognized mesh format.").raise();
      if (fmt == MeshFormat::MFEM)
      {
         getObject() = Rodin::Mesh<Traits::Serial>(mfem::Mesh(is, 0, 1, getFixOrientation()));
         return {true, {}};
      }
      else
      {
         return {false, IO::Error{"Cannot determine MFEM mesh format version."}};
      }
   }

   IO::Status MeshLoader<MeshFormat::GMSH, Traits::Serial>::load(std::istream& is)
   {
      assert(is);
      auto fmt = getMeshFormat(is);
      if (!fmt.has_value())
         Alert::Exception("Unrecognized mesh format.").raise();
      if (fmt == MeshFormat::GMSH)
      {
         getObject() = Rodin::Mesh<Traits::Serial>(mfem::Mesh(is, 0, 1, getFixOrientation()));
         return {true, {}};
      }
      else
      {
         return {false, IO::Error{"Cannot determine Gmsh mesh format version."}};
      }
   }

   IO::Status MeshLoader<MeshFormat::MEDIT, Traits::Serial>::load(std::istream& is)
   {
      int spaceDim = 0;
      auto fmt = getMeshFormat(is, false);
      if (!fmt.has_value())
         Alert::Exception("Unrecognized mesh format.").raise();
      if (fmt == MeshFormat::MEDIT)
      {
         mfem::Mesh mfemMesh;

         std::string line;

         std::map<Medit::EntityKeyword, std::optional<std::istream::pos_type>> pos;
         std::map<Medit::EntityKeyword, size_t> count;

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
            else if (Medit::EntityKeywordMap.left.count(kw))
            {
               auto ent = Medit::EntityKeywordMap.left.at(kw);
               std::getline(is, line);
               pos[ent] = is.tellg();
               std::istringstream lss(line);
               lss >> count[ent];
            }
         }

         // Infer type of mesh
         bool isSurfaceMesh;
         if (spaceDim == 3 && (
                  count.count(Medit::EntityKeyword::Tetrahedra) == 0
                  || count.at(Medit::EntityKeyword::Tetrahedra) == 0))
         {
            // It's a surface mesh
            mfemMesh = mfem::Mesh(
                  spaceDim - 1, // Dimension
                  count.at(Medit::EntityKeyword::Vertices),
                  count.at(Medit::EntityKeyword::Triangles),
                  count.at(Medit::EntityKeyword::Edges),
                  spaceDim // Space dimension
                  );
            isSurfaceMesh = true;
         }
         else if (spaceDim == 2 && (
                  count.count(Medit::EntityKeyword::Triangles) == 0
                  || count.at(Medit::EntityKeyword::Triangles) == 0))
         {
            // It's a surface mesh
            mfemMesh = mfem::Mesh(
                  spaceDim - 1, // Dimension
                  count.at(Medit::EntityKeyword::Vertices),
                  count.at(Medit::EntityKeyword::Edges),
                  spaceDim // Space dimension
                  );
            isSurfaceMesh = true;
         }
         else if (spaceDim == 3)
         {
            // It's a volume mesh
            mfemMesh = mfem::Mesh(
                  spaceDim,
                  count.at(Medit::EntityKeyword::Vertices),
                  count.at(Medit::EntityKeyword::Tetrahedra),
                  count.at(Medit::EntityKeyword::Triangles),
                  spaceDim
                  );
            isSurfaceMesh = false;
         }
         else if (spaceDim == 2)
         {
            // It's a volume mesh
            mfemMesh = mfem::Mesh(
                  spaceDim,
                  count.at(Medit::EntityKeyword::Vertices),
                  count.at(Medit::EntityKeyword::Triangles),
                  count.at(Medit::EntityKeyword::Edges),
                  spaceDim
                  );
            isSurfaceMesh = false;
         }
         else
         {
            Alert::Exception("Unhandled case.").raise();
         }

         is.clear();
         for (const auto& [kw, ent] : Medit::EntityKeywordMap)
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
                        case Medit::EntityKeyword::Vertices:
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
                              mfemMesh.AddVertex(coords[0], coords[1], coords[2]);
                              // We ignore the reference
                           }
                           break;
                        }
                        case Medit::EntityKeyword::Triangles:
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
                                 Alert::Exception("Unhandled case.").raise();
                              }
                           }
                           break;
                        }
                        case Medit::EntityKeyword::Tetrahedra:
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
                              mfemMesh.AddTet(v1 - 1, v2 - 1, v3 - 1, v4 - 1, ref);
                           }
                           break;
                        }
                        case Medit::EntityKeyword::Edges:
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

         mfemMesh.FinalizeMesh(0, getFixOrientation());

         getObject() = Rodin::Mesh<Traits::Serial>(std::move(mfemMesh));

         return {true, {}};
      }
      else
      {
         return {false, IO::Error{"Cannot determine Medit mesh format version."}};
      }
   }
}
