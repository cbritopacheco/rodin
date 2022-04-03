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
#include "ScalarSolution3D.h"

namespace Rodin::External::MMG
{
  // ---- ScalarSolution3D --------------------------------------------------
  ScalarSolution3D::ScalarSolution3D(MMG5_pSol sol, Mesh3D& mesh)
    : ScalarSolution(sol),
      m_mesh(mesh)
  {}

  ScalarSolution3D::ScalarSolution3D(Mesh3D& mesh)
    : ScalarSolution(),
      m_mesh(mesh)
  {
    getHandle()->ver  = 2;
    getHandle()->size = 1;
    getHandle()->type = MMG5_Scalar;
    getHandle()->dim  = mesh.getDimension(); // Supported on 3D mesh

    getHandle()->np  = mesh.count(Mesh3D::Vertex);
    getHandle()->npi = getHandle()->npi;
    getHandle()->npmax = std::max(static_cast<int>(1.5 * getHandle()->np), MMG3D_NPMAX);
    MMG5_SAFE_CALLOC(
          getHandle()->m, (getHandle()->size * (getHandle()->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
  }

  ScalarSolution3D::ScalarSolution3D(const ScalarSolution3D& other)
    : ScalarSolution(other),
      m_mesh(other.m_mesh)
  {}

  ScalarSolution3D::ScalarSolution3D(ScalarSolution3D&& other)
    : ScalarSolution(std::move(other)),
      m_mesh(other.m_mesh)
  {}

  ScalarSolution3D&
  ScalarSolution3D::operator=(const ScalarSolution3D& other)
  {
     if (this != &other)
     {
        m_mesh = other.m_mesh;
        MMG5_Sol_Swap(getHandle(), ScalarSolution3D(other).getHandle());
     }

     return *this;
  }

  IncompleteScalarSolution3D ScalarSolution3D::load(
        const boost::filesystem::path& filename)
  {
     IncompleteScalarSolution3D res;
     Load_MMG5_Sol(filename, 2, res.getHandle());
     return res;
  }

  void ScalarSolution3D::save(const boost::filesystem::path& filename)
  {
     if (!getHandle()->np || !getHandle()->m)
     {
        Alert::Exception(
              "Failed to write ScalarSolution3D to file. No data!").raise();
     }

     if (!MMG3D_saveSol(m_mesh.get().getHandle(), getHandle(), filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  ScalarSolution3D& ScalarSolution3D::setMesh(Mesh3D& mesh)
  {
     m_mesh = mesh;
     return *this;
  }

  const Mesh3D& ScalarSolution3D::getMesh() const
  {
     return m_mesh;
  }

  Mesh3D& ScalarSolution3D::getMesh()
  {
     return m_mesh;
  }

  // ---- IncompleteScalarSolution3D -------------------------------------------
  IncompleteScalarSolution3D::IncompleteScalarSolution3D()
    : IncompleteSolutionBase()
  {
    getHandle()->ver  = 2;
    getHandle()->size = 1;
    getHandle()->type = MMG5_Scalar;
    getHandle()->dim  = 3;
  }
}
