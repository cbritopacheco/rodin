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
#include <mmg2d/mmg2d.h>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Warning.h"

#include "Utility.h"

#include "Mesh2D.h"

namespace Rodin::External::MMG
{
  Mesh2D Mesh2D::load(const std::filesystem::path& filename)
  {
    Mesh2D mesh;
    if (!MMG2D_loadMesh(mesh.getHandle(), filename.c_str()))
    {
       Alert::Exception(
             "Failed to open file for reading: " + filename.string()).raise();
    }
    return mesh;
  }

  void Mesh2D::save(const std::filesystem::path& filename)
  {
    if (!MMG2D_saveMesh(getHandle(), filename.c_str()))
    {
       Alert::Exception(
             "Failed to open file for writing: " + filename.string()).raise();
    }
  }

  Mesh2D::Mesh2D()
  {
    m_mesh = NULL;

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

    m_mesh->dim   = 2;
    m_mesh->ver   = 2;
    m_mesh->nsols = 0;

    MMG2D_Set_commonFunc();
    MMG2D_Init_parameters(m_mesh);
    MMG2D_Init_fileNames(m_mesh, NULL);

    // Verbosity is high when in Debug
    MMG2D_Set_iparameter(
          getHandle(),
          nullptr,
          MMG2D_IPARAM_verbose,
          VERBOSITY_LEVEL);
  }

  Mesh2D::Mesh2D(const Mesh2D& other)
     : Mesh2D()
  {
    MMG5_Mesh_Copy(other.getHandle(), getHandle());
  }

  Mesh2D::~Mesh2D()
  {
    if (m_mesh)
      MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &m_mesh, MMG5_ARG_end);
  }

  Mesh2D::Mesh2D(Mesh2D&& other)
  {
    m_mesh = other.m_mesh;
    other.m_mesh = nullptr;
  }

  int Mesh2D::count(Mesh2D::Entity e) const
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
        Alert::Exception("Unknown Mesh2D::Entity").raise();
    }
    return 0;
  }

  MMG5_pMesh& Mesh2D::getHandle()
  {
    return m_mesh;
  }

  const MMG5_pMesh& Mesh2D::getHandle() const
  {
    return m_mesh;
  }
}
