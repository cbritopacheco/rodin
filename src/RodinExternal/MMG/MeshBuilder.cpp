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
    m_requiredEdgeIndex = std::move(other.m_requiredEdgeIndex);
    m_requiredVertexIndex = std::move(other.m_requiredVertexIndex);
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

  Mesh::Builder& Mesh::Builder::requiredEdge(Index edgeIdx)
  {
    m_requiredEdgeIndex.insert(edgeIdx);
    return *this;
  }

  Mesh::Builder& Mesh::Builder::requiredVertex(Index vertexIdx)
  {
    m_requiredVertexIndex.insert(vertexIdx);
    return *this;
  }

  MMG::Mesh Mesh::Builder::finalize()
  {
    MMG::Mesh res;
    res.Parent::operator=(Parent::Builder::finalize());
    res.m_cornerIndex = std::move(m_cornerIndex);
    res.m_ridgeIndex = std::move(m_ridgeIndex);
    res.m_requiredEdgeIndex = std::move(m_requiredEdgeIndex);
    res.m_requiredVertexIndex = std::move(m_requiredVertexIndex);
    return res;
  }
}

