/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert/Exception.h"

#include "Mesh.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"

namespace Rodin::External::MMG
{
  Mesh& Mesh::setCorner(Index vertexIdx)
  {
    m_cornerIndex.insert(vertexIdx);
    return *this;
  }

  Mesh& Mesh::setRidge(Index edgeIdx)
  {
    m_ridgeIndex.insert(edgeIdx);
    return *this;
  }

  void Mesh::save(
     const boost::filesystem::path& filename, IO::FileFormat fmt, size_t precison) const
  {
    if (fmt == IO::FileFormat::MEDIT)
    {
      std::ofstream os(filename.c_str());
      if (!os)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to open " << filename << " for writing."
          << Alert::Raise;
      }
      os.precision(precison);
      MMG::MeshPrinter printer(*this);
      printer.print(os);
    }
    else
    {
      Parent::save(filename, fmt, precison);
    }
  }

  Mesh& Mesh::load(
     const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    if (fmt == IO::FileFormat::MEDIT)
    {
      std::ifstream in(filename.c_str());
      if (!in)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to open " << filename << " for reading."
          << Alert::Raise;
      }
      MMG::MeshLoader loader(*this);
      loader.load(in);
    }
    else
    {
      Parent::load(filename, fmt);
    }
    return *this;
  }
}

