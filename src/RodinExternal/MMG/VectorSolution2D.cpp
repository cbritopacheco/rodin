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
#include "VectorSolution2D.h"

namespace Rodin::External::MMG
{
  // ---- VectorSolution2D --------------------------------------------------
  VectorSolution2D::VectorSolution2D(MMG5_pSol sol, Mesh2D& mesh)
    : VectorSolution(sol),
      m_mesh(mesh)
  {}

  VectorSolution2D::VectorSolution2D(Mesh2D& mesh)
    : VectorSolution(),
      m_mesh(mesh)
  {
    getHandle()->ver  = 2;
    getHandle()->size = 2;
    getHandle()->type = MMG5_Vector;
    getHandle()->dim  = mesh.getDimension(); // Supported on 2D mesh

    getHandle()->np  = mesh.count(Mesh2D::Vertex);
    getHandle()->npi = getHandle()->npi;
    getHandle()->npmax = std::max(static_cast<int>(1.5 * getHandle()->np), MMG2D_NPMAX);
    MMG5_SAFE_CALLOC(
          getHandle()->m, (getHandle()->size * (getHandle()->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
  }

  VectorSolution2D::VectorSolution2D(const VectorSolution2D& other)
    : VectorSolution(other),
      m_mesh(other.m_mesh)
  {}

  VectorSolution2D::VectorSolution2D(VectorSolution2D&& other)
    : VectorSolution(std::move(other)),
      m_mesh(other.m_mesh)
  {}

  VectorSolution2D&
  VectorSolution2D::operator=(const VectorSolution2D& other)
  {
     if (this != &other)
     {
        m_mesh = other.m_mesh;
        MMG5_Sol_Swap(getHandle(), VectorSolution2D(other).getHandle());
     }

     return *this;
  }

  IncompleteVectorSolution2D VectorSolution2D::load(
        const boost::filesystem::path& filename)
  {
     IncompleteVectorSolution2D res;
     Load_MMG5_Sol(filename, 2, res.getHandle());
     return res;
  }

  void VectorSolution2D::save(const boost::filesystem::path& filename)
  {
     if (!getHandle()->np || !getHandle()->m)
     {
        Alert::Exception(
              "Failed to write VectorSolution2D to file. No data!").raise();
     }

     if (!MMG2D_saveSol(m_mesh.get().getHandle(), getHandle(), filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  VectorSolution2D& VectorSolution2D::setMesh(Mesh2D& mesh)
  {
     m_mesh = mesh;
     return *this;
  }

  const Mesh2D& VectorSolution2D::getMesh() const
  {
     return m_mesh;
  }

  Mesh2D& VectorSolution2D::getMesh()
  {
     return m_mesh;
  }

  // ---- IncompleteVectorSolution2D -------------------------------------------
  IncompleteVectorSolution2D::IncompleteVectorSolution2D()
    : IncompleteSolutionBase()
  {
    getHandle()->ver  = 2;
    getHandle()->size = 2;
    getHandle()->type = MMG5_Vector;
    getHandle()->dim  = 2;
  }
}
