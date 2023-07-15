/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mesh.h"

namespace Rodin::External::MMG
{
  Mesh::Builder& Mesh::Builder::operator=(Mesh::Builder&& other)
  {
    Parent::Builder::operator=(std::move(other));
    m_cornerIndex = std::move(other.m_cornerIndex);
    m_ridgeIndex = std::move(other.m_ridgeIndex);
    return *this;
  }

  Mesh::Builder& Mesh::Builder::corner(Index vertexIdx)
  {
    m_cornerIndex.insert(vertexIdx);
    return *this;
  }

  Mesh::Builder& Mesh::Builder::ridge(Index edgeIdx)
  {
    m_ridgeIndex.insert(edgeIdx);
    return *this;
  }

  Mesh Mesh::Builder::finalize()
  {
    Mesh res;
    res.Parent::operator=(Parent::Builder::finalize());
    res.m_cornerIndex = std::move(m_cornerIndex);
    res.m_ridgeIndex = std::move(m_ridgeIndex);
    return res;
  }
}

