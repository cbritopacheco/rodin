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
    // Copy all non-pointer fields
    *m_mesh = *other.m_mesh;

    // Fields not used in two dimensions
    m_mesh->adjt = nullptr;
    m_mesh->adjapr = nullptr;
    m_mesh->xpoint = nullptr;
    m_mesh->tetra = nullptr;
    m_mesh->ipar = nullptr;
    m_mesh->xtetra = nullptr;
    m_mesh->prism = nullptr;
    m_mesh->xprism = nullptr;

    // Pointer fields
    if (other.m_mesh->adja)
    {
      MMG5_SAFE_CALLOC(m_mesh->adja, 4 * other.m_mesh->nemax + 5, int, /* No op */);
      std::copy(
          other.m_mesh->adja, other.m_mesh->adja + 4 * other.m_mesh->nemax + 5,
          m_mesh->adja);
    }
    else
      m_mesh->adja = nullptr;

    if (other.m_mesh->adjq)
    {
      MMG5_SAFE_CALLOC(m_mesh->adjq, 4 * other.m_mesh->nquad + 5, int, /* No op */);
      std::copy(
          other.m_mesh->adjq, other.m_mesh->adjq + 4 * other.m_mesh->nquad + 5,
          m_mesh->adjq);
    }
    else
      m_mesh->adjq = nullptr;

    // Why MMG fields are 1-indexed is beyond me, but it does not pose so much
    // difficulty. We will copy from the 0-index, if it contains garbage then
    // it should not matter because MMG code should only access memory from the
    // 1-index.
    if (other.m_mesh->point)
    {
      MMG5_SAFE_CALLOC(m_mesh->point, other.m_mesh->npmax + 1,
          MMG5_Point, /* No op */);
      std::copy(
          other.m_mesh->point, other.m_mesh->point + other.m_mesh->npmax + 1,
          m_mesh->point);
    }
    else
      m_mesh->point = nullptr;

    if (other.m_mesh->tria)
    {
      MMG5_SAFE_CALLOC(m_mesh->tria, other.m_mesh->ntmax + 1, MMG5_Tria, /* No op */);
      std::copy(
          other.m_mesh->tria, other.m_mesh->tria + other.m_mesh->ntmax + 1,
          m_mesh->tria);
    }
    else
      m_mesh->tria = nullptr;

    if (other.m_mesh->quadra)
    {
      MMG5_SAFE_CALLOC(m_mesh->quadra, other.m_mesh->nquad + 1, MMG5_Quad, /* No op */);
      std::copy(
          other.m_mesh->quadra, other.m_mesh->quadra + other.m_mesh->nquad + 1,
          m_mesh->quadra);
    }
    else
      m_mesh->quadra = nullptr;

    if (other.m_mesh->edge)
    {
      MMG5_SAFE_CALLOC(m_mesh->edge, other.m_mesh->namax + 1, MMG5_Edge, /* No op */);
      std::copy(
          other.m_mesh->edge, other.m_mesh->edge + other.m_mesh->namax + 1,
          m_mesh->edge);
    }
    else
      m_mesh->edge = nullptr;

    m_mesh->info = other.m_mesh->info;
    m_mesh->info.br = nullptr;
    if (other.m_mesh->info.mat)
    {
      MMG5_SAFE_CALLOC(m_mesh->info.mat, other.m_mesh->info.nmat + 1,
          MMG5_Mat, /* No op */);
      std::copy(
          other.m_mesh->info.mat, other.m_mesh->info.mat + other.m_mesh->info.nmat + 1,
          m_mesh->info.mat);
    }
    else
      m_mesh->info.mat = nullptr;

    if (other.m_mesh->info.invmat.lookup)
    {
      MMG5_SAFE_CALLOC(m_mesh->info.invmat.lookup,
          other.m_mesh->info.invmat.size, int, /* No op */);
      std::copy(
          other.m_mesh->info.invmat.lookup,
          other.m_mesh->info.invmat.lookup + other.m_mesh->info.invmat.size,
          m_mesh->info.invmat.lookup
          );
    }
    else
      m_mesh->info.invmat.lookup = nullptr;

    MMG5_SAFE_CALLOC(m_mesh->namein, std::strlen(other.m_mesh->namein),
        char, /* No op */);
    std::strcpy(m_mesh->namein, other.m_mesh->namein);
    MMG5_SAFE_CALLOC(m_mesh->nameout, std::strlen(other.m_mesh->nameout),
        char, /* No op */);
    std::strcpy(m_mesh->nameout, other.m_mesh->nameout);
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
