/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "SubMesh.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Local>::Builder&
  SubMesh<Context::Local>::Builder::initialize(const Mesh<Context>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();
    m_parent = parent;
    m_build.initialize(sdim);
    m_s2ps.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  SubMesh<Context::Local>::Builder&
  SubMesh<Context::Local>::Builder::include(size_t d, Index parentIdx)
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
      const auto [it, inserted] = m_s2ps[0].right.insert(std::pair<Index, Index>{ parentVertex, childVertex });
      if (inserted) // Vertex was not already in the map
      {
        m_s2ps[0].left.push_back(parentVertex);
        childPolytope.coeffRef(i) = childVertex;
        m_sidx[0] += 1;
      }
      else // Vertex was already in the map
      {
        childPolytope.coeffRef(i) = it->second;
      }
    }
    // Add polytope with original geometry and new vertex ordering
    build.polytope(conn.getGeometry(d, parentIdx), childPolytope);
    const Index childIdx = m_sidx[d];
    const auto [it, inserted] = m_s2ps[d].right.insert(std::pair<Index, Index>{ parentIdx, childIdx });
    // Add polytope information
    if (inserted) // Polytope was already in the map
    {
      m_s2ps[d].left.push_back(parentIdx);
      build.attribute(std::pair<Index, Index>{ d, childIdx }, parent.getAttribute(d, parentIdx));
      m_sidx[d] += 1;
    }
    else
    {
      build.attribute({ d, it->second }, parent.getAttribute(d, parentIdx));
    }
    m_dimension = std::max(m_dimension, d);
    return *this;
  }

  SubMesh<Context::Local>::Builder&
  SubMesh<Context::Local>::Builder::include(size_t d, const IndexSet& indices)
  {
    for (const Index parentIdx : indices)
      include(d, parentIdx);
    return *this;
  }

  SubMesh<Context::Local> SubMesh<Context::Local>::Builder::finalize()
  {
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    const size_t nodes = m_sidx[0];

    // Build the mesh object.
    m_build.nodes(nodes);
    for (const Index& pIdx : m_s2ps[0].left)
      m_build.vertex(parent.getVertexCoordinates(pIdx));

    // Build the connectivity for the submesh from the parent mesh.
    auto& conn = m_build.getConnectivity();
    for (size_t d = 0; d < m_s2ps.size(); d++)
    {
      for (size_t dp = 0; dp < m_s2ps.size(); dp++)
      {
        if (d == m_dimension && dp == 0)
          continue;
        const auto& pInc = parent.getConnectivity().getIncidence(d, dp);
        if (pInc.size() > 0)
        {
          assert(m_s2ps[d].left.size() == m_s2ps[d].right.size());
          Incidence cInc(m_s2ps[d].left.size());
          for (size_t i = 0; i < m_s2ps[d].left.size(); ++i)
          {
            const Index cIdx = i;
            const Index pIdx = m_s2ps[d].left[i];
            cInc[cIdx].reserve(pInc[pIdx].size());
            for (const Index p : pInc[pIdx])
            {
              auto find = m_s2ps[dp].right.find(p);
              if (find != m_s2ps[dp].right.end())
                cInc[cIdx].insert_unique(find->second);
            }
          }
          // Manually set the incidence
          conn.setIncidence({ d, dp }, std::move(cInc));
        }
      }
    }
    // Finalize construction
    SubMesh res(parent);
    res.Parent::operator=(m_build.finalize());
    res.m_s2ps = std::move(m_s2ps);
    return res;
  }
}
