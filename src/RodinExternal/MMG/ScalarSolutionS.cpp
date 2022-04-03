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
#include "ScalarSolutionS.h"

namespace Rodin::External::MMG
{
  // ---- ScalarSolutionS --------------------------------------------------
  ScalarSolutionS::ScalarSolutionS(MMG5_pSol sol, MeshS& mesh)
    : ScalarSolution(sol),
      m_mesh(mesh)
  {}

  ScalarSolutionS::ScalarSolutionS(MeshS& mesh)
    : ScalarSolution(),
      m_mesh(mesh)
  {
    getHandle()->ver  = 2;
    getHandle()->size = 1;
    getHandle()->type = MMG5_Scalar;
    getHandle()->dim  = mesh.getDimension(); // Supported on S mesh

    getHandle()->np  = mesh.count(MeshS::Vertex);
    getHandle()->npi = getHandle()->npi;
    getHandle()->npmax = std::max(static_cast<int>(1.5 * getHandle()->np), MMGS_NPMAX);
    MMG5_SAFE_CALLOC(
          getHandle()->m, (getHandle()->size * (getHandle()->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
  }

  ScalarSolutionS::ScalarSolutionS(const ScalarSolutionS& other)
    : ScalarSolution(other),
      m_mesh(other.m_mesh)
  {}

  ScalarSolutionS::ScalarSolutionS(ScalarSolutionS&& other)
    : ScalarSolution(std::move(other)),
      m_mesh(other.m_mesh)
  {}

  ScalarSolutionS&
  ScalarSolutionS::operator=(const ScalarSolutionS& other)
  {
     if (this != &other)
     {
        m_mesh = other.m_mesh;
        MMG5_Sol_Swap(getHandle(), ScalarSolutionS(other).getHandle());
     }

     return *this;
  }

  IncompleteScalarSolutionS ScalarSolutionS::load(
        const boost::filesystem::path& filename)
  {
     IncompleteScalarSolutionS res;
     Load_MMG5_Sol(filename, 3, res.getHandle());
     return res;
  }

  void ScalarSolutionS::save(const boost::filesystem::path& filename)
  {
     if (!getHandle()->np || !getHandle()->m)
     {
        Alert::Exception(
              "Failed to write ScalarSolutionS to file. No data!").raise();
     }

     if (!MMGS_saveSol(m_mesh.get().getHandle(), getHandle(), filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  ScalarSolutionS& ScalarSolutionS::setMesh(MeshS& mesh)
  {
     m_mesh = mesh;
     return *this;
  }

  const MeshS& ScalarSolutionS::getMesh() const
  {
     return m_mesh;
  }

  MeshS& ScalarSolutionS::getMesh()
  {
     return m_mesh;
  }

  // ---- IncompleteScalarSolutionS -------------------------------------------
  IncompleteScalarSolutionS::IncompleteScalarSolutionS()
    : IncompleteSolutionBase()
  {
    getHandle()->ver  = 2;
    getHandle()->size = 1;
    getHandle()->type = MMG5_Scalar;
    getHandle()->dim  = 3;
  }
}
