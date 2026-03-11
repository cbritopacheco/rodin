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
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    const auto& conn   = parent.getConnectivity();

    // Fast path: entity already included -> only ensure attribute (optional), then return.
    if (auto it = m_s2ps[d].right.find(parentIdx); it != m_s2ps[d].right.end())
    {
      if (auto attr = parent.getAttribute(d, parentIdx))
        m_build.attribute({ d, it->second }, *attr);
      return *this;
    }

    // First time we see this entity: create child index now
    const Index childIdx = m_sidx[d];
    m_s2ps[d].left.push_back(parentIdx);
    m_s2ps[d].right.insert(std::pair<Index, Index>{ parentIdx, childIdx });
    ++m_sidx[d];

    // Remap vertices and emit polytope ONCE
    const auto& parentPolytope = conn.getPolytope(d, parentIdx);
    IndexArray childPolytope(parentPolytope.size());

    for (size_t i = 0; i < static_cast<size_t>(childPolytope.size()); ++i)
    {
      const Index pv = parentPolytope(i);

      // Vertex map check (same idea can be optimized further; see section 2)
      auto [vit, vinserted] = m_s2ps[0].right.insert(std::pair<Index, Index>{ pv, m_sidx[0] });
      if (vinserted)
      {
        m_s2ps[0].left.push_back(pv);
        ++m_sidx[0];
      }
      childPolytope.coeffRef(i) = vit->second;
    }

    m_build.polytope(conn.getGeometry(d, parentIdx), std::move(childPolytope));

    if (auto attr = parent.getAttribute(d, parentIdx))
      m_build.attribute({ d, childIdx }, *attr);

    m_dimension = std::max(m_dimension, d);
    return *this;
  }

  SubMesh<Context::Local>::Builder&
  SubMesh<Context::Local>::Builder::include(size_t d, const IndexVector& indices)
  {
    for (const Index parentIdx : indices)
      this->include(d, parentIdx);
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
                cInc[cIdx].push_back(find->second);
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
