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
    return m_s2ds[d];
  }

  size_t Shard::Builder::getPolytopeCount(size_t d) const
  {
    assert(m_s2ds[d].left.size() == m_s2ds[d].right.size());
    return m_s2ds[d].left.size();
  }

  Shard::Builder& Shard::Builder::initialize(const Mesh<Context>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();
    m_parent = parent;
    m_build.initialize(sdim);
    m_s2ds.resize(dim + 1);
    m_flags.resize(dim + 1);
    m_owner.resize(dim + 1);
    m_halo.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  UnorderedMap<Index, Index>& Shard::Builder::getOwner(size_t d)
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  UnorderedMap<Index, IndexSet>& Shard::Builder::getHalo(size_t d)
  {
    assert(d < m_halo.size());
    return m_halo[d];
  }

  const UnorderedMap<Index, Index>& Shard::Builder::getOwner(size_t d) const
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  const UnorderedMap<Index, IndexSet>& Shard::Builder::getHalo(size_t d) const
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
      const auto [it, inserted] = m_s2ds[0].right.insert(std::pair<Index, Index>{ parentIdx, childIdx });
      if (inserted) // Vertex was not already in the map
      {
        m_s2ds[0].left.push_back(parentIdx);
        auto attr = parent.getAttribute(0, parentIdx);
        if (attr)
          build.attribute({ 0, childIdx }, *attr);
        m_flags[0].push_back(flags);
        m_sidx[0] += 1;
      }
      res.first = it->second;
      res.second = inserted;
    }
    else
    {
      const Index childIdx = m_sidx[d];
      const auto [it, inserted] = m_s2ds[d].right.insert(std::pair<Index, Index>{ parentIdx, childIdx });
      assert(parentPolytope.size());
      IndexArray childPolytope(parentPolytope.size());
      // Add polytope information
      if (inserted) // Polytope was not already in the map
      {
        m_s2ds[d].left.push_back(parentIdx);
        auto attr = parent.getAttribute(d, parentIdx);
        if (attr)
          build.attribute({ d, childIdx }, *attr);
        m_flags[d].push_back(flags);
        m_sidx[d] += 1;
        for (size_t i = 0; i < static_cast<size_t>(childPolytope.size()); i++)
        {
          const Index parentVertex = parentPolytope(i);
          const auto find = m_s2ds[0].right.find(parentVertex);
          if (find != m_s2ds[0].right.end()) // Vertex is in the map
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
      res.first = it->second;
      res.second = inserted;
    }
    m_dimension = std::max(m_dimension, d);
    return res;
  }

  Shard Shard::Builder::finalize()
  {
    assert(m_parent.has_value());
    const auto& parent = m_parent->get();
    const auto& pconn  = parent.getConnectivity();
    const size_t D     = parent.getDimension();

    auto& cconn = m_build.getConnectivity();

    // --------------------------------------------------------------------------
    // Rebuild all incidences on the local shard by filtering parent incidences
    // through the shard-local polytope maps.
    // --------------------------------------------------------------------------
    for (size_t d = 0; d <= D; ++d)
    {
      const auto& dmap = m_s2ds[d];
      const size_t nd  = dmap.left.size();
      if (nd == 0)
        continue;

      for (size_t dp = 0; dp <= D; ++dp)
      {
        const auto& pInc = pconn.getIncidence(d, dp);
        if (pInc.empty())
          continue;

        const auto& dpmap = m_s2ds[dp].right;
        Incidence cInc(nd);

        for (Index ci = 0; ci < static_cast<Index>(nd); ++ci)
        {
          const Index pi = dmap.left[ci];
          const auto& pinc = pInc.at(pi);

          auto& cinc = cInc[ci];
          cinc.reserve(pinc.size());

          for (const Index pj : pinc)
          {
            auto it = dpmap.find(pj);
            if (it != dpmap.end())
              cinc.push_back(it->second);
          }
        }

        cconn.setIncidence({ d, dp }, std::move(cInc));
      }
    }

    // --------------------------------------------------------------------------
    // Rebuild vertex coordinates in shard-local order.
    // Vertices are emitted in the same order as m_s2ds[0].left.
    // --------------------------------------------------------------------------
    const auto& vmap = m_s2ds[0].left;
    m_build.nodes(vmap.size());
    for (const Index pvid : vmap)
      m_build.vertex(parent.getVertexCoordinates(pvid));

    // --------------------------------------------------------------------------
    // Finalize the local mesh and move shard metadata.
    // --------------------------------------------------------------------------
    Shard res;
    res.Parent::operator=(m_build.finalize());
    res.m_s2ds  = std::move(m_s2ds);
    res.m_flags = std::move(m_flags);
    res.m_owner = std::move(m_owner);
    res.m_halo  = std::move(m_halo);

    return res;
  }

  Shard::Shard(const Shard& other)
    : Parent(other),
      m_s2ds(other.m_s2ds),
      m_flags(other.m_flags),
      m_owner(other.m_owner),
      m_halo(other.m_halo)
  {}

  Shard::Shard(Shard&& other)
    : Parent(std::move(other)),
      m_s2ds(std::move(other.m_s2ds)),
      m_flags(std::move(other.m_flags)),
      m_owner(std::move(other.m_owner)),
      m_halo(std::move(other.m_halo))
  {}

  Shard& Shard::operator=(Shard&& other)
  {
    Parent::operator=(std::move(other));
    m_s2ds = std::move(other.m_s2ds);
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

  const UnorderedMap<Index, Index>& Shard::getOwner(size_t d) const
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  const UnorderedMap<Index, IndexSet>& Shard::getHalo(size_t d) const
  {
    assert(d < m_halo.size());
    return m_halo[d];
  }

  const Shard::PolytopeMap& Shard::getPolytopeMap(size_t d) const
  {
    return m_s2ds[d];
  }

  UnorderedMap<Index, Index>& Shard::getOwner(size_t d)
  {
    assert(d < m_owner.size());
    return m_owner[d];
  }

  UnorderedMap<Index, IndexSet>& Shard::getHalo(size_t d)
  {
    assert(d < m_halo.size());
    return m_halo[d];
  }

  Shard::PolytopeMap& Shard::getPolytopeMap(size_t d)
  {
    return m_s2ds[d];
  }

  const std::vector<Shard::Flags>& Shard::getFlags(size_t d) const
  {
    return m_flags[d];
  }

  std::vector<Shard::Flags>& Shard::getFlags(size_t d)
  {
    return m_flags[d];
  }
}

