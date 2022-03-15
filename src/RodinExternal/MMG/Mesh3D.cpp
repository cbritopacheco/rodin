/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <limits>
#include <cstring>
#include <fstream>

#include <libmmgcommon.h>
#include <mmg3d/mmg3d.h>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Warning.h"

#include "Utility.h"

#include "Mesh3D.h"

namespace Rodin::External::MMG
{
   Mesh3D::Mesh3D()
   {
       m_mesh = nullptr;

       auto calloc =
          [this]()
          {
             // Returns false on fail
             MMG5_SAFE_CALLOC(m_mesh, 1, MMG5_Mesh, return false);
             // Returns true on success
             return true;
          };

       if (!calloc())
          Alert::Exception("Failed to allocate memory for the mesh").raise();

       m_mesh->dim   = 3;
       m_mesh->ver   = 2;
       m_mesh->nsols = 0;

       MMG3D_Set_commonFunc();
       MMG3D_Init_parameters(m_mesh);
       MMG3D_Init_fileNames(m_mesh, nullptr);

       // Verbosity is high when in Debug
       MMG3D_Set_iparameter(
             getHandle(),
             nullptr,
             MMG3D_IPARAM_verbose,
             VERBOSITY_LEVEL);
   }

   Mesh3D::Mesh3D(const Mesh3D& other)
      : Mesh3D()
   {
      MMG5_Mesh_Copy(other.getHandle(), getHandle());
   }


   Mesh3D::Mesh3D(Mesh3D&& other)
   {
       m_mesh = other.m_mesh;
       other.m_mesh = nullptr;
   }

   Mesh3D::~Mesh3D()
   {
       if (m_mesh)
         MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &m_mesh, MMG5_ARG_end);
   }

   Mesh3D Mesh3D::load(const std::filesystem::path& filename)
   {
     Mesh3D mesh;
     if (!MMG3D_loadMesh(mesh.getHandle(), filename.c_str()))
     {
        Alert::Exception(
              "Failed to open file for reading: " + filename.string()).raise();
     }
     return mesh;
   }

   void Mesh3D::save(const std::filesystem::path& filename)
   {
     if (!MMG3D_saveMesh(getHandle(), filename.c_str()))
     {
        Alert::Exception(
              "Failed to open file for writing: " + filename.string()).raise();
     }
   }

  int Mesh3D::count(Mesh3D::Entity e) const
  {
    switch (e)
    {
      case Entity::Vertex:
        return getHandle()->np;
      case Entity::Edge:
        return getHandle()->na;
      case Entity::Triangle:
        return getHandle()->nt;
      case Entity::Tetrahedra:
        return getHandle()->ne;
      default:
        Alert::Exception("Unknown Mesh3D::Entity").raise();
    }
    return 0;
  }

   MMG5_pMesh& Mesh3D::getHandle()
   {
      return m_mesh;
   }

   const MMG5_pMesh& Mesh3D::getHandle() const
   {
      return m_mesh;
   }
}
