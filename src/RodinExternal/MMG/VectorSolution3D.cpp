/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <string>
#include <cstdio>
#include <cstring>

#include "Utility.h"
#include "VectorSolution3D.h"

namespace Rodin::External::MMG
{
  // ---- VectorSolution3D --------------------------------------------------
  VectorSolution3D::VectorSolution3D(MMG5_pSol sol, Mesh3D& mesh)
    : VectorSolution(sol),
      m_mesh(mesh)
  {}

  VectorSolution3D::VectorSolution3D(Mesh3D& mesh)
    : VectorSolution(),
      m_mesh(mesh)
  {
    getHandle()->ver  = 2;
    getHandle()->size = 3;
    getHandle()->type = MMG5_Vector;
    getHandle()->dim  = mesh.getDimension(); // Supported on 3D mesh

    getHandle()->np  = mesh.count(Mesh3D::Vertex);
    getHandle()->npi = getHandle()->npi;
    getHandle()->npmax = std::max(static_cast<int>(1.5 * getHandle()->np), MMG3D_NPMAX);
    MMG5_SAFE_CALLOC(
          getHandle()->m, (getHandle()->size * (getHandle()->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
  }

  VectorSolution3D::VectorSolution3D(const VectorSolution3D& other)
    : VectorSolution(other),
      m_mesh(other.m_mesh)
  {}

  VectorSolution3D::VectorSolution3D(VectorSolution3D&& other)
    : VectorSolution(std::move(other)),
      m_mesh(other.m_mesh)
  {}

  VectorSolution3D&
  VectorSolution3D::operator=(const VectorSolution3D& other)
  {
     if (this != &other)
     {
        m_mesh = other.m_mesh;
        MMG5_Sol_Swap(getHandle(), VectorSolution3D(other).getHandle());
     }

     return *this;
  }

  VectorSolution3D& VectorSolution3D::load(const boost::filesystem::path& filename)
  {
     Load_MMG5_Sol(filename, 3, getHandle());
     return *this;
  }

  void VectorSolution3D::save(const boost::filesystem::path& filename)
  {
     if (!getHandle()->np || !getHandle()->m)
     {
        Alert::Exception(
              "Failed to write VectorSolution3D to file. No data!").raise();
     }

     if (!MMG3D_saveSol(m_mesh.get().getHandle(), getHandle(), filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  VectorSolution3D& VectorSolution3D::setMesh(Mesh3D& mesh)
  {
     m_mesh = mesh;
     return *this;
  }

  const Mesh3D& VectorSolution3D::getMesh() const
  {
     return m_mesh;
  }

  Mesh3D& VectorSolution3D::getMesh()
  {
     return m_mesh;
  }
}

