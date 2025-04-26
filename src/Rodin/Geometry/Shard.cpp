#include "Shard.h"

namespace Rodin::Geometry
{
  Shard::Builder& Shard::Builder::initialize(const Mesh<Context>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();
    m_parent = parent;
    m_build.initialize(sdim);
    m_s2ps.resize(dim + 1);
    m_ghosts.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  Shard::Builder& Shard::Builder::include(size_t d, Index parentIdx)
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
    if (inserted) // Polytope was not already in the map
    {
      build.attribute({ d, childIdx }, parent.getAttribute(d, parentIdx));
      m_sidx[d] += 1;
    }
    else
    {
      build.attribute({ d, it->get_left() }, parent.getAttribute(d, parentIdx));
    }
    m_dimension = std::max(m_dimension, d);
    return *this;
  }

  Shard::Builder& Shard::Builder::ghost(size_t d, Index parentIdx)
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
        m_ghosts[0].insert(childVertex);
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
    if (inserted) // Polytope was not already in the map
    {
      build.attribute({ d, childIdx }, parent.getAttribute(d, parentIdx));
      m_sidx[d] += 1;
      m_ghosts[d].insert(childIdx);
    }
    else
    {
      build.attribute({ d, it->get_left() }, parent.getAttribute(d, parentIdx));
    }
    m_dimension = std::max(m_dimension, d);
    return *this;
  }

  Shard::Builder& Shard::Builder::include(size_t d, const IndexSet& indices)
  {
    for (const Index parentIdx : indices)
      include(d, parentIdx);
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
      if (!m_ghosts[cellDim].contains(cIdx))
        originals.push_back(it->get_right());
    }

    const auto& nbrs = pconn.getIncidence(cellDim, cellDim);
    for (Index p : originals)
    {
      for (Index nb : nbrs.at(p))
      {
        if (m_s2ps[cellDim].right.find(nb) == m_s2ps[cellDim].right.end())
          ghost(cellDim, nb);
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
            ghost(d - 1, sub);
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
    res.m_s2ps   = std::move(m_s2ps);
    res.m_ghosts = std::move(m_ghosts);
    return res;
  }

  Shard::Shard(const Shard& other)
    : Parent(other),
      m_s2ps(other.m_s2ps),
      m_ghosts(other.m_ghosts)
  {}

  Shard::Shard(Shard&& other)
    : Parent(std::move(other)),
      m_s2ps(std::move(other.m_s2ps)),
      m_ghosts(std::move(other.m_ghosts))
  {}

  Shard& Shard::operator=(Shard&& other)
  {
    Parent::operator=(std::move(other));
    m_s2ps = std::move(other.m_s2ps);
    m_ghosts = std::move(other.m_ghosts);
    return *this;
  }

  bool Shard::isGhost(size_t d, Index idx) const
  {
    return m_ghosts[d].contains(idx);
  }

  const boost::bimap<Index, Index>& Shard::getPolytopeMap(size_t d) const
  {
    return m_s2ps[d];
  }

  Real Shard::getVolume() const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(3); !it.end(); ++it)
    {
      if (!isGhost(3, it->getIndex()))
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real Shard::getVolume(Attribute attr) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(3); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (!isGhost(3, it->getIndex()))
          totalVolume += it->getMeasure();
      }
    }
    return totalVolume;
  }

  Real Shard::getVolume(const FlatSet<Attribute>& attrs) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(3); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
      {
        if (!isGhost(3, it->getIndex()))
          totalVolume += it->getMeasure();
      }
    }
    return totalVolume;
  }

  Real Shard::getPerimeter() const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (!isGhost(it->getDimension(), it->getIndex()))
        totalPerimeter += it->getMeasure();
    }
    return totalPerimeter;
  }

  Real Shard::getPerimeter(Attribute attr) const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (!isGhost(it->getDimension(), it->getIndex()))
          totalPerimeter += it->getMeasure();
      }
    }
    return totalPerimeter;
  }

  Real Shard::getPerimeter(const FlatSet<Attribute>& attrs) const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
      {
        if (!isGhost(it->getDimension(), it->getIndex()))
          totalPerimeter += it->getMeasure();
      }
    }
    return totalPerimeter;
  }

  Real Shard::getArea() const
  {
    Real totalArea = 0;
    for (auto it = getPolytope(2); !it.end(); ++it)
    {
      if (!isGhost(2, it->getIndex()))
        totalArea += it->getMeasure();
    }
    return totalArea;
  }

  Real Shard::getArea(Attribute attr) const
  {
    Real totalArea = 0;
    for (auto it = getPolytope(2); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (!isGhost(2, it->getIndex()))
          totalArea += it->getMeasure();
      }
    }
    return totalArea;
  }

  Real Shard::getArea(const FlatSet<Attribute>& attrs) const
  {
    Real totalArea = 0;
    for (auto it = getPolytope(2); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
      {
        if (!isGhost(2, it->getIndex()))
          totalArea += it->getMeasure();
      }
    }
    return totalArea;
  }

  Real Shard::getMeasure(size_t d) const
  {
    Real totalMeasure = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (!isGhost(d, it->getIndex()))
        totalMeasure += it->getMeasure();
    }
    return totalMeasure;
  }

  Real Shard::getMeasure(size_t d, Attribute attr) const
  {
    Real totalMeasure = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (!isGhost(d, it->getIndex()))
          totalMeasure += it->getMeasure();
      }
    }
    return totalMeasure;
  }

  Real Shard::getMeasure(size_t d, const FlatSet<Attribute>& attrs) const
  {
    Real totalMeasure = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
      {
        if (!isGhost(d, it->getIndex()))
          totalMeasure += it->getMeasure();
      }
    }
    return totalMeasure;
  }
}

