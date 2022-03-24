#include <libmmgcommon.h>
#include <mmgs/mmgs.h>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Warning.h"

#include "Utility.h"

#include "MeshS.h"

namespace Rodin::External::MMG
{
   MeshS::MeshS()
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

       MMGS_Set_commonFunc();
       MMGS_Init_parameters(m_mesh);
       MMGS_Init_fileNames(m_mesh, nullptr);

       // Verbosity is high when in Debug
       MMGS_Set_iparameter(
             getHandle(),
             nullptr,
             MMGS_IPARAM_verbose,
             VERBOSITY_LEVEL);
   }

   MeshS::MeshS(const MeshS& other)
      : MeshS()
   {
      MMG5_Mesh_Copy(other.getHandle(), getHandle());
   }


   MeshS::MeshS(MeshS&& other)
   {
       m_mesh = other.m_mesh;
       other.m_mesh = nullptr;
   }

   MeshS::~MeshS()
   {
       if (m_mesh)
         MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &m_mesh, MMG5_ARG_end);
   }

   MeshS MeshS::load(const std::filesystem::path& filename)
   {
     MeshS mesh;
     if (!MMGS_loadMesh(mesh.getHandle(), filename.c_str()))
     {
        Alert::Exception(
              "Failed to open file for reading: " + filename.string()).raise();
     }
     return mesh;
   }

   void MeshS::save(const std::filesystem::path& filename)
   {
     if (!MMGS_saveMesh(getHandle(), filename.c_str()))
     {
        Alert::Exception(
              "Failed to open file for writing: " + filename.string()).raise();
     }
   }

  int MeshS::count(MeshS::Entity e) const
  {
    switch (e)
    {
      case Entity::Vertex:
        return getHandle()->np;
      case Entity::Edge:
        return getHandle()->na;
      case Entity::Triangle:
        return getHandle()->nt;
      default:
        Alert::Exception("Unknown MeshS::Entity").raise();
    }
    return 0;
  }

   MMG5_pMesh& MeshS::getHandle()
   {
      return m_mesh;
   }

   const MMG5_pMesh& MeshS::getHandle() const
   {
      return m_mesh;
   }
}
