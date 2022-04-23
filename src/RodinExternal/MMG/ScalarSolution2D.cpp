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
#include "ScalarSolution2D.h"

namespace Rodin::External::MMG
{
  // ---- ScalarSolution2D --------------------------------------------------
  ScalarSolution2D::ScalarSolution2D(MMG5_pSol sol, Mesh2D& mesh)
    : ScalarSolution(sol),
      m_mesh(mesh)
  {}

  ScalarSolution2D::ScalarSolution2D(Mesh2D& mesh)
    : ScalarSolution(),
      m_mesh(mesh)
  {
    getHandle()->ver  = 2;
    getHandle()->size = 1;
    getHandle()->type = MMG5_Scalar;
    getHandle()->dim  = mesh.getDimension(); // Supported on 2D mesh

    getHandle()->np  = mesh.count(Mesh2D::Vertex);
    getHandle()->npi = getHandle()->npi;
    getHandle()->npmax = std::max(static_cast<int>(1.5 * getHandle()->np), MMG2D_NPMAX);
    MMG5_SAFE_CALLOC(
          getHandle()->m, (getHandle()->size * (getHandle()->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
  }

  ScalarSolution2D::ScalarSolution2D(const ScalarSolution2D& other)
    : ScalarSolution(other),
      m_mesh(other.m_mesh)
  {}

  ScalarSolution2D::ScalarSolution2D(ScalarSolution2D&& other)
    : ScalarSolution(std::move(other)),
      m_mesh(other.m_mesh)
  {}

  ScalarSolution2D&
  ScalarSolution2D::operator=(const ScalarSolution2D& other)
  {
     if (this != &other)
     {
        m_mesh = other.m_mesh;
        MMG5_Sol_Swap(getHandle(), ScalarSolution2D(other).getHandle());
     }

     return *this;
  }

  ScalarSolution2D& ScalarSolution2D::load(const boost::filesystem::path& filename)
  {
     Load_MMG5_Sol(filename, 2, getHandle());
     return *this;
  }

  void ScalarSolution2D::save(const boost::filesystem::path& filename)
  {
     if (!getHandle()->np || !getHandle()->m)
     {
        Alert::Exception(
              "Failed to write ScalarSolution2D to file. No data!").raise();
     }

     if (!MMG2D_saveSol(m_mesh.get().getHandle(), getHandle(), filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  ScalarSolution2D& ScalarSolution2D::setMesh(Mesh2D& mesh)
  {
     m_mesh = mesh;
     return *this;
  }

  const Mesh2D& ScalarSolution2D::getMesh() const
  {
     return m_mesh;
  }

  Mesh2D& ScalarSolution2D::getMesh()
  {
     return m_mesh;
  }
}
