#include "Shard.h"
#include "Rodin/Alert/Identifier.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"

namespace Rodin::Geometry
{
  const Shard::PolytopeMap& Shard::Builder::getPolytopeMap(size_t d) const
  {
    return m_s2ds[d];
  }

  size_t Shard::Builder::getPolytopeCount(size_t d) const
  {
    assert(m_s2ds[d].left.size() == m_s2ds[d].right.size());
    return m_s2ds[d].left.size();
  }

  Shard::Builder::Builder()
    : m_dimension(0),
      m_sdim(0),
      m_mode(Mode::None)
  {}

  Shard::Builder& Shard::Builder::initialize(const Mesh<Context>& parent)
  {
    const size_t dim = parent.getDimension();
    const size_t sdim = parent.getSpaceDimension();

    m_mode = Mode::Parent;

    m_parent = parent;
    m_build = Mesh<Context>::Builder();
    m_build.initialize(sdim);

    m_s2ds.clear();
    m_state.clear();
    m_owner.clear();
    m_halo.clear();
    m_sidx.clear();
    m_vertices.clear();

    m_s2ds.resize(dim + 1);
    m_state.resize(dim + 1);
    m_owner.resize(dim + 1);
    m_halo.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);

    m_vertices.setDimension(static_cast<std::uint8_t>(sdim));

    m_dimension = dim;
    return *this;
  }

  Shard::Builder& Shard::Builder::initialize(size_t dimension, size_t sdim)
  {
    m_mode = Mode::Direct;
    m_parent.reset();

    m_dimension = dimension;
    m_sdim = sdim;

    m_build = Mesh<Context>::Builder();
    m_build.initialize(sdim);

    m_s2ds.clear();
    m_state.clear();
    m_owner.clear();
    m_halo.clear();
    m_sidx.clear();

    m_s2ds.resize(dimension + 1);
    m_state.resize(dimension + 1);
    m_owner.resize(dimension + 1);
    m_halo.resize(dimension + 1);
    m_sidx.resize(dimension + 1, 0);

    m_connectivity.initialize(sdim);
    m_attributes.initialize(sdim);
    m_transformations.initialize(sdim);

    return *this;
  }

  Index Shard::Builder::vertex(
      Index globalIdx,
      const Math::SpatialPoint& x,
      const State& state)
  {
    if (m_mode != Mode::Direct)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "vertex(...) is only valid in direct-construction mode."
        << Alert::Raise;
    }

    auto& vmap = m_s2ds[0];
    if (auto it = vmap.right.find(globalIdx); it != vmap.right.end())
      return it->second;

    const Index localIdx = m_sidx[0]++;
    const auto [it, inserted] = vmap.right.emplace(globalIdx, localIdx);
    (void) it;
    assert(inserted);

    vmap.left.push_back(globalIdx);
    m_state[0].push_back(state);

    m_vertices.push_back(x);
    m_attributes.resize(0, m_sidx[0]);
    m_transformations.resize(0, m_sidx[0]);

    return localIdx;
  }

  Index Shard::Builder::polytope(
      size_t d,
      Index globalIdx,
      Polytope::Type g,
      const IndexArray& vs,
      const State& state)
  {
    if (m_mode != Mode::Direct)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "This method is only valid in direct-construction mode."
        << Alert::Raise;
    }

    if (d == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Use " << Alert::Identifier::Function("vertex") << " to insert zero-dimensional polytopes."
        << Alert::Raise;
    }

    assert(d < m_s2ds.size());

    auto& map = m_s2ds[d];
    if (auto it = map.right.find(globalIdx); it != map.right.end())
      return it->second;

    for (size_t i = 0; i < static_cast<size_t>(vs.size()); ++i)
    {
      if (vs[i] >= m_sidx[0])
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Vertex local index " << vs[i] << " is out of bounds."
          << Alert::Raise;
      }
    }

    const Index localIdx = m_sidx[d]++;
    const auto [it, inserted] = map.right.emplace(globalIdx, localIdx);
    (void)it;
    assert(inserted);

    map.left.push_back(globalIdx);
    m_state[d].push_back(state);

    m_connectivity.polytope(g, vs);
    m_attributes.resize(d, m_sidx[d]);
    m_transformations.resize(d, m_sidx[d]);

    return localIdx;
  }

  Shard::Builder& Shard::Builder::setOwner(size_t d, Index localIdx, Index ownerRank)
  {
    assert(d < m_owner.size());
    assert(localIdx < m_s2ds[d].left.size());
    m_owner[d][localIdx] = ownerRank;
    return *this;
  }

  Shard::Builder& Shard::Builder::halo(size_t d, Index localIdx, Index neighborRank)
  {
    assert(d < m_halo.size());
    assert(localIdx < m_s2ds[d].left.size());
    m_halo[d][localIdx].insert(neighborRank);
    return *this;
  }

  Shard::Builder& Shard::Builder::attribute(
      const std::pair<size_t, Index>& p,
      const Optional<Attribute>& attr)
  {
    const auto& [d, i] = p;
    assert(d < m_s2ds.size());
    assert(i < m_s2ds[d].left.size());

    if (m_mode == Mode::Parent)
    {
      m_build.attribute(p, attr);
    }
    else
    {
      m_attributes.set(p, attr);
    }

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
      const std::pair<size_t, Index>& p, const State& state)
  {
    if (m_mode != Mode::Parent || !m_parent.has_value())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "This method requires parent-based initialization."
        << Alert::Raise;
    }

    const auto& [d, parentIdx] = p;
    assert(m_parent.has_value());
    const auto& parent = m_parent.value().get();
    const auto& conn = parent.getConnectivity();
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
        m_state[0].push_back(state);
        m_sidx[0] += 1;
      }
      res.first = it->second;
      res.second = inserted;
    }
    else
    {
      const auto& parentPolytope = conn.getPolytope(d, parentIdx);
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
        m_state[d].push_back(state);
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
    if (m_mode == Mode::Parent)
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
    }
    else
    {
      m_vertices.setDimension(static_cast<std::uint8_t>(m_sdim));
      m_connectivity.nodes(m_s2ds[0].left.size());
      m_build.setVertices(std::move(m_vertices))
             .setConnectivity(std::move(m_connectivity))
             .setAttributeIndex(std::move(m_attributes))
             .setTransformationIndex(std::move(m_transformations));
    }

    // --------------------------------------------------------------------------
    // Finalize the local mesh and move shard metadata.
    // --------------------------------------------------------------------------
    Shard res;
    res.Parent::operator=(m_build.finalize());
    res.m_s2ds  = std::move(m_s2ds);
    res.m_state = std::move(m_state);
    res.m_owner = std::move(m_owner);
    res.m_halo  = std::move(m_halo);

    return res;
  }

  Shard::Shard(const Shard& other)
    : Parent(other),
      m_s2ds(other.m_s2ds),
      m_state(other.m_state),
      m_owner(other.m_owner),
      m_halo(other.m_halo)
  {}

  Shard::Shard(Shard&& other)
    : Parent(std::move(other)),
      m_s2ds(std::move(other.m_s2ds)),
      m_state(std::move(other.m_state)),
      m_owner(std::move(other.m_owner)),
      m_halo(std::move(other.m_halo))
  {}

  Shard& Shard::operator=(Shard&& other)
  {
    Parent::operator=(std::move(other));
    m_s2ds = std::move(other.m_s2ds);
    m_state = std::move(other.m_state);
    m_owner = std::move(other.m_owner);
    m_halo = std::move(other.m_halo);
    return *this;
  }

  bool Shard::isOwned(size_t d, Index idx) const
  {
    return m_state[d][idx] == State::Owned;
  }

  bool Shard::isShared(size_t d, Index idx) const
  {
    return m_state[d][idx] == State::Shared;
  }

  bool Shard::isGhost(size_t d, Index idx) const
  {
    return m_state[d][idx] == State::Ghost;
  }

  bool Shard::isLocal(size_t d, Index idx) const
  {
    const auto s = m_state[d][idx];
    return s == State::Owned || s == State::Shared;
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

  const std::vector<Shard::State>& Shard::getState(size_t d) const
  {
    return m_state[d];
  }

  std::vector<Shard::State>& Shard::getState(size_t d)
  {
    return m_state[d];
  }
}

