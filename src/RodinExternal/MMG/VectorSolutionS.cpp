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
#include "VectorSolutionS.h"

namespace Rodin::External::MMG
{
  // ---- VectorSolutionS --------------------------------------------------
  VectorSolutionS::VectorSolutionS(MMG5_pSol sol, MeshS& mesh)
    : VectorSolution(sol),
      m_mesh(mesh)
  {}

  VectorSolutionS::VectorSolutionS(MeshS& mesh)
    : VectorSolution(),
      m_mesh(mesh)
  {
    getHandle()->ver  = 2;
    getHandle()->size = 3;
    getHandle()->type = MMG5_Vector;
    getHandle()->dim  = mesh.getDimension(); // Supported on S mesh

    getHandle()->np  = mesh.count(MeshS::Vertex);
    getHandle()->npi = getHandle()->npi;
    getHandle()->npmax = std::max(static_cast<int>(1.5 * getHandle()->np), MMGS_NPMAX);
    MMG5_SAFE_CALLOC(
          getHandle()->m, (getHandle()->size * (getHandle()->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
  }

  VectorSolutionS::VectorSolutionS(const VectorSolutionS& other)
    : VectorSolution(other),
      m_mesh(other.m_mesh)
  {}

  VectorSolutionS::VectorSolutionS(VectorSolutionS&& other)
    : VectorSolution(std::move(other)),
      m_mesh(other.m_mesh)
  {}

  VectorSolutionS&
  VectorSolutionS::operator=(const VectorSolutionS& other)
  {
     if (this != &other)
     {
        m_mesh = other.m_mesh;
        MMG5_Sol_Swap(getHandle(), VectorSolutionS(other).getHandle());
     }

     return *this;
  }

  VectorSolutionS& VectorSolutionS::load(const boost::filesystem::path& filename)
  {
     Load_MMG5_Sol(filename, 3, getHandle());
     return *this;
  }

  void VectorSolutionS::save(const boost::filesystem::path& filename)
  {
     if (!getHandle()->np || !getHandle()->m)
     {
        Alert::Exception(
          "Failed to write VectorSolutionS to file. No data!").raise();
     }

     if (!MMGS_saveSol(m_mesh.get().getHandle(), getHandle(), filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  VectorSolutionS& VectorSolutionS::setMesh(MeshS& mesh)
  {
     m_mesh = mesh;
     return *this;
  }

  const MeshS& VectorSolutionS::getMesh() const
  {
     return m_mesh;
  }

  MeshS& VectorSolutionS::getMesh()
  {
     return m_mesh;
  }
}
