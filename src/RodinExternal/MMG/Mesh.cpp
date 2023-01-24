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
  Mesh& Mesh::corner(int vertexIdx)
  {
    m_corners.insert(vertexIdx);
    return *this;
  }

  Mesh& Mesh::ridge(int edgeIdx)
  {
    assert(edgeIdx >= 0);
    m_ridges.insert(edgeIdx);
    return *this;
  }

  Mesh& Mesh::edge(const std::pair<int, int>& endpoints, int ref)
  {
    m_edges.push_back(Edge{endpoints, ref});
    return *this;
  }

  void Mesh::save(
     const boost::filesystem::path& filename, IO::FileFormat fmt, size_t precison) const
  {
    if (fmt == IO::FileFormat::MEDIT)
    {
      std::ofstream os(filename.c_str());
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

