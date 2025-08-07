#include "Shard.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"

namespace Rodin::Geometry
{
  Shard::Flags const Shard::Flags::None(0b00);

  Shard::Flags const Shard::Flags::Owned(0b01);

  Shard::Flags const Shard::Flags::Ghost(0b10);

  const Shard::PolytopeMap& Shard::Builder::getPolytopeMap(size_t d) const
  {
    return m_s2ps[d];
  }

  size_t Shard::Builder::getPolytopeCount(size_t d) const
  {
    assert(m_s2ps[d].left.size() == m_s2ps[d].right.size());
    return m_s2ps[d].left.size();
  }

  Shard::Builder& Shard::Builder::initialize(const Mesh<Context>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();
    m_parent = parent;
    m_build.initialize(sdim);
    m_s2ps.resize(dim + 1);
    m_flags.resize(dim + 1);
    m_owner.resize(dim + 1);
    m_halo.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  FlatMap<Index, Index>& Shard::Builder::getOwner(size_t d)
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  FlatMap<Index, IndexSet>& Shard::Builder::getHalo(size_t d)
  {
    assert(d < m_halo.size());
    return m_halo[d];
  }

  const FlatMap<Index, Index>& Shard::Builder::getOwner(size_t d) const
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  const FlatMap<Index, IndexSet>& Shard::Builder::getHalo(size_t d) const
  {
    assert(d < m_halo.size());
    return m_halo[d];
  }

  std::pair<Index, Boolean> Shard::Builder::include(
      const std::pair<size_t, Index>& p, const Flags& flags)
  {
    const auto& [d, parentIdx] = p;
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    const auto& conn = parent.getConnectivity();
    const auto& parentPolytope = conn.getPolytope(d, parentIdx);
    auto& build = m_build;

    std::pair<Index, Boolean> res;

    if (d == 0)
    {
      const Index childIdx = m_sidx[0];
      const auto [it, inserted] = m_s2ps[0].right.insert({ parentIdx, childIdx });
      if (inserted) // Vertex was not already in the map
      {
        m_s2ps[0].left.push_back(parentIdx);
        build.attribute({ 0, childIdx }, parent.getAttribute(0, parentIdx));
        m_flags[0].push_back(flags);
        m_sidx[0] += 1;
      }
      res = { it->second, inserted };
    }
    else
    {
      const Index childIdx = m_sidx[d];
      const auto [it, inserted] = m_s2ps[d].right.insert({ parentIdx, childIdx });
      assert(parentPolytope.size());
      IndexArray childPolytope(parentPolytope.size());
      // Add polytope information
      if (inserted) // Polytope was not already in the map
      {
        m_s2ps[d].left.push_back(parentIdx);
        build.attribute({ d, childIdx }, parent.getAttribute(d, parentIdx));
        m_flags[d].push_back(flags);
        m_sidx[d] += 1;
        for (size_t i = 0; i < static_cast<size_t>(childPolytope.size()); i++)
        {
          const Index parentVertex = parentPolytope.coeff(i);
          const auto find = m_s2ps[0].right.find(parentVertex);
          if (find != m_s2ps[0].right.end()) // Vertex is in the map
          {
            childPolytope.coeffRef(i) = find->second;
          }
          else // Vertex is not in the map
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Vertex " << parentVertex << " of polytope " << parentIdx
              << " of dimension " << d << " is not in the map."
              << Alert::Raise;
          }
        }
        // Add polytope with original geometry and new vertex ordering
        build.polytope(conn.getGeometry(d, parentIdx), childPolytope);
      }
      res = { it->second, inserted };
    }
    m_dimension = std::max(m_dimension, d);
    return res;
  }

  Shard Shard::Builder::finalize()
  {
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    auto& conn = m_build.getConnectivity();
    const size_t cellDim = parent.getDimension();
    const auto& pconn = parent.getConnectivity();
    for (size_t d = 0; d <= cellDim; ++d)
    {
      for (size_t dp = 0; dp <= cellDim; ++dp)
      {
        const auto& pInc = pconn.getIncidence(d, dp);
        if (pInc.empty())
          continue;
        Incidence cInc(m_s2ps[d].left.size());
        for (size_t i = 0; i < m_s2ps[d].left.size(); ++i)
        {
          const Index& cIdx = i;
          const Index& pIdx = m_s2ps[d].left[i];
          cInc[cIdx].reserve(pInc.at(pIdx).size());
          for (Index pNbr : pInc.at(pIdx))
          {
            const auto found = m_s2ps[dp].right.find(pNbr);
            if (found != m_s2ps[dp].right.end())
              cInc[cIdx].insert_unique(found->second);
          }
        }
        conn.setIncidence({ d, dp }, std::move(cInc));
      }
    }

    m_build.nodes(m_sidx[0]);
    for (const Index& pIdx: m_s2ps[0].left)
      m_build.vertex(parent.getVertexCoordinates(pIdx));

    Shard res;
    res.Parent::operator=(m_build.finalize());
    res.m_s2ps = std::move(m_s2ps);
    res.m_flags = std::move(m_flags);
    res.m_owner = std::move(m_owner);
    res.m_halo = std::move(m_halo);
    return res;
  }

  Shard::Shard(const Shard& other)
    : Parent(other),
      m_s2ps(other.m_s2ps),
      m_flags(other.m_flags),
      m_owner(other.m_owner),
      m_halo(other.m_halo)
  {}

  Shard::Shard(Shard&& other)
    : Parent(std::move(other)),
      m_s2ps(std::move(other.m_s2ps)),
      m_flags(std::move(other.m_flags)),
      m_owner(std::move(other.m_owner)),
      m_halo(std::move(other.m_halo))
  {}

  Shard& Shard::operator=(Shard&& other)
  {
    Parent::operator=(std::move(other));
    m_s2ps = std::move(other.m_s2ps);
    m_flags = std::move(other.m_flags);
    m_owner = std::move(other.m_owner);
    m_halo = std::move(other.m_halo);
    return *this;
  }

  bool Shard::isGhost(size_t d, Index idx) const
  {
    return m_flags[d][idx] & Shard::Flags::Ghost;
  }

  bool Shard::isOwned(size_t d, Index idx) const
  {
    return m_flags[d][idx] & Shard::Flags::Owned;
  }

  const FlatMap<Index, Index>& Shard::getOwner(size_t d) const
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  const FlatMap<Index, IndexSet>& Shard::getHalo(size_t d) const
  {
    assert(d < m_halo.size());
    return m_halo[d];
  }

  const Shard::PolytopeMap& Shard::getPolytopeMap(size_t d) const
  {
    return m_s2ps[d];
  }
}

