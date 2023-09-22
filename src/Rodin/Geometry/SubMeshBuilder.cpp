/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "SubMesh.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::initialize(const Mesh<Context::Serial>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();
    m_parent = parent;
    m_build.initialize(sdim);
    m_s2ps.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::include(size_t d, Index parentIdx)
  {
    auto& build = m_build;
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    const auto& conn = parent.getConnectivity();
    const auto& parentPolytope = conn.getPolytope(d, parentIdx);
    IndexArray childPolytope(parentPolytope.size());
    assert(childPolytope.size() >= 0);
    for (size_t i = 0; i < static_cast<size_t>(childPolytope.size()); i++)
    {
      const Index parentVertex = parentPolytope.coeff(i);
      const Index childVertex = m_sidx[0];
      const auto [it, inserted] = m_s2ps[0].left.insert({ childVertex, parentVertex });
      if (inserted) // Vertex was not already in the map
      {
        childPolytope.coeffRef(i) = childVertex;
        m_sidx[0] += 1;
      }
      else // Vertex was already in the map
      {
        childPolytope.coeffRef(i) = it->get_left();
      }
    }
    // Add polytope with original geometry and new vertex ordering
    build.polytope(conn.getGeometry(d, parentIdx), childPolytope);
    const Index childIdx = m_sidx[d];
    const auto [it, inserted] = m_s2ps[d].left.insert({ childIdx, parentIdx });
    // Add polytope information
    if (inserted) // Polytope was already in the map
    {
      build.attribute({ d, childIdx }, parent.getAttribute(d, parentIdx));
      m_sidx[d] += 1;
    }
    else
    {
      build.attribute({ d, it->get_left() }, parent.getAttribute(d, parentIdx));
    }

    return *this;
  }

  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::include(size_t d, const IndexSet& indices)
  {
    for (const Index parentIdx : indices)
      include(d, parentIdx);
    return *this;
  }

  SubMesh<Context::Serial> SubMesh<Context::Serial>::Builder::finalize()
  {
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    const size_t nodes = m_sidx[0];
    m_build.nodes(nodes);
    for (auto it = m_s2ps[0].left.begin(); it != m_s2ps[0].left.end(); ++it)
      m_build.vertex(parent.getVertexCoordinates(it->get_right()));
    SubMesh res(parent);
    res.Parent::operator=(m_build.finalize());
    res.m_s2ps = std::move(m_s2ps);
    return res;
  }
}
