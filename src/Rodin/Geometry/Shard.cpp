#include "Shard.h"

namespace Rodin::Geometry
{
  Shard::Flags const Shard::Flags::None(0b00);

  Shard::Flags const Shard::Flags::Owned(0b01);

  Shard::Flags const Shard::Flags::Ghost(0b10);

  Shard::Builder& Shard::Builder::initialize(const Mesh<Context>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();
    m_parent = parent;
    m_build.initialize(sdim);
    m_s2ps.resize(dim + 1);
    m_flags.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  Shard::Builder& Shard::Builder::include(size_t d, Index parentIdx, const Flags& flags)
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
      const auto [it, inserted] = m_s2ps[0].right.insert({ parentVertex, childVertex });
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
    const auto [it, inserted] = m_s2ps[d].right.insert({ parentIdx, childIdx });
    // Add polytope information
    if (inserted) // Polytope was not already in the map
    {
      build.attribute({ d, childIdx }, parent.getAttribute(d, parentIdx));
      m_flags[d].insert({ childIdx, flags });
      m_sidx[d] += 1;
    }
    else
    {
      build.attribute({ d, it->get_left() }, parent.getAttribute(d, parentIdx));
      m_flags[d][it->get_left()] |= flags;
    }
    m_dimension = std::max(m_dimension, d);
    return *this;
  }

  Shard Shard::Builder::finalize()
  {
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    auto&       conn   = m_build.getConnectivity();
    const size_t cellDim = parent.getDimension();
    const auto& pconn   = parent.getConnectivity();

    std::vector<Index> originals;
    for (auto it = m_s2ps[cellDim].left.begin(); it != m_s2ps[cellDim].left.end(); ++it)
    {
      Index cIdx = it->get_left();
      const auto& flags = m_flags[cellDim].at(cIdx);
      if (flags & Shard::Flags::Ghost)
      {
        include(cellDim, it->get_right(), Shard::Flags::Ghost);
        continue;
      }
    }

    const auto& nbrs = pconn.getIncidence(cellDim, cellDim);
    for (Index p : originals)
    {
      for (Index nb : nbrs.at(p))
      {
        if (m_s2ps[cellDim].right.find(nb) == m_s2ps[cellDim].right.end())
          include(cellDim, nb, Shard::Flags::Ghost);
      }
    }

    for (size_t d = cellDim; d > 0; --d)
    {
      const auto& down = pconn.getIncidence(d, d - 1);
      if (down.size() == 0) continue;
      for (auto it = m_s2ps[d].left.begin(); it != m_s2ps[d].left.end(); ++it)
      {
        Index pParent = it->get_right();
        for (Index sub : down.at(pParent))
          if (m_s2ps[d - 1].right.find(sub) == m_s2ps[d - 1].right.end())
            include(d - 1, sub, Shard::Flags::Ghost);
      }
    }

    for (size_t d = 0; d <= cellDim; ++d)
    {
      for (size_t dp = 0; dp <= cellDim; ++dp)
      {
        const auto& pInc = pconn.getIncidence(d, dp);
        if (pInc.empty()) continue;

        Incidence cInc(m_s2ps[d].size());
        for (auto it = m_s2ps[d].left.begin(); it != m_s2ps[d].left.end(); ++it)
        {
          Index cIdx = it->get_left();
          Index pIdx = it->get_right();
          cInc[cIdx].reserve(pInc.at(pIdx).size());
          for (Index pNbr : pInc.at(pIdx))
          {
            auto found = m_s2ps[dp].right.find(pNbr);
            if (found != m_s2ps[dp].right.end())
              cInc[cIdx].insert_unique(found->get_left());
          }
        }
        conn.setIncidence({ d, dp }, std::move(cInc));
      }
    }

    m_build.nodes(m_sidx[0]);
    for (auto it = m_s2ps[0].left.begin(); it != m_s2ps[0].left.end(); ++it)
      m_build.vertex(parent.getVertexCoordinates(it->get_right()));

    Shard res;
    res.Parent::operator=(m_build.finalize());
    res.m_s2ps = std::move(m_s2ps);
    res.m_flags = std::move(m_flags);
    return res;
  }

  Shard::Shard(const Shard& other)
    : Parent(other),
      m_s2ps(other.m_s2ps),
      m_flags(other.m_flags)
  {}

  Shard::Shard(Shard&& other)
    : Parent(std::move(other)),
      m_s2ps(std::move(other.m_s2ps)),
      m_flags(std::move(other.m_flags))
  {}

  Shard& Shard::operator=(Shard&& other)
  {
    Parent::operator=(std::move(other));
    m_s2ps = std::move(other.m_s2ps);
    m_flags = std::move(other.m_flags);
    return *this;
  }

  bool Shard::isGhost(size_t d, Index idx) const
  {
    return m_flags[d].at(idx) & Shard::Flags::Ghost;
  }

  const Shard::PolytopeMap& Shard::getPolytopeMap(size_t d) const
  {
    return m_s2ps[d];
  }
}

