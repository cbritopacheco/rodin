#include <boost/dynamic_bitset.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeTransformation.h"
#include "Rodin/Types.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  MPIMesh::Builder::Builder(const Context::MPI& context)
    : m_context(context)
  {}

  MPIMesh::Builder& MPIMesh::Builder::initialize(Shard&& shard)
  {
    m_shard = std::move(shard);
    return *this;
  }

  MPIMesh MPIMesh::Builder::finalize()
  {
    MPIMesh mesh(m_context);
    mesh.m_shard = std::move(m_shard);
    return mesh;
  }

  MPIMesh& MPIMesh::scale(Real c)
  {
    this->getShard().scale(c);
    return *this;
  }

  void MPIMesh::flush()
  {
    this->getShard().flush();
  }

  bool MPIMesh::isSubMesh() const
  {
    return false;
  }

  size_t MPIMesh::getDimension() const
  {
    return this->getShard().getDimension();
  }

  size_t MPIMesh::getSpaceDimension() const
  {
    return this->getShard().getSpaceDimension();
  }

  const Context::MPI& MPIMesh::getContext() const
  {
    return m_context;
  }

  size_t MPIMesh::getPolytopeCount(size_t d) const
  {
    size_t localCount = 0;
    const auto& shard = getShard();
    const auto& ctx = getContext();
    const boost::mpi::communicator& comm = ctx.getCommunicator();
    for (size_t i = 0; i < shard.getPolytopeCount(d); ++i)
    {
      if (shard.isOwned(d, i))
        ++localCount;
    }
    return boost::mpi::all_reduce(comm, localCount, std::plus<size_t>());
  }

  size_t MPIMesh::getPolytopeCount(Polytope::Type g) const
  {
    const size_t d = Polytope::Traits(g).getDimension();
    size_t localCount = 0;
    const auto& shard = getShard();
    const auto& ctx = getContext();
    const boost::mpi::communicator& comm = ctx.getCommunicator();
    for (size_t i = 0; i < shard.getPolytopeCount(d); ++i)
    {
      if (shard.getGeometry(d, i) == g)
      {
        if (shard.isOwned(d, i))
          ++localCount;
      }
    }
    return boost::mpi::all_reduce(comm, localCount, std::plus<size_t>());
  }

  Shard& MPIMesh::getShard()
  {
    return m_shard;
  }

  const Shard& MPIMesh::getShard() const
  {
    return m_shard;
  }

  Optional<Index> MPIMesh::getLocalIndex(size_t dimension, Index globalIdx) const
  {
    const auto& map = getShard().getPolytopeMap(dimension).right;
    auto it = map.find(globalIdx);
    if (it == map.end())
      return std::nullopt;
    return it->second;
  }

  Index MPIMesh::getGlobalIndex(size_t dimension, Index localIdx) const
  {
    const auto& shard = getShard();
    const auto& pm = shard.getPolytopeMap(dimension);
    return pm.left.at(localIdx);
  }

  CellIterator MPIMesh::getCell() const
  {
    return this->getCell(0);
  }

  FaceIterator MPIMesh::getFace() const
  {
    return this->getFace(0);
  }

  VertexIterator MPIMesh::getVertex() const
  {
    return this->getVertex(0);
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension) const
  {
    return this->getPolytope(dimension, 0);
  }

  CellIterator MPIMesh::getCell(Index localIdx) const
  {
    return CellIterator(*this, BoundedIndexGenerator(localIdx, this->getShard().getCellCount()));
  }

  FaceIterator MPIMesh::getFace(Index localIdx) const
  {
    return FaceIterator(*this, BoundedIndexGenerator(localIdx, this->getShard().getFaceCount()));
  }

  VertexIterator MPIMesh::getVertex(Index localIdx) const
  {
    return VertexIterator(*this, BoundedIndexGenerator(localIdx, this->getShard().getVertexCount()));
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index localIdx) const
  {
    return PolytopeIterator(dimension, *this, BoundedIndexGenerator(localIdx, this->getShard().getPolytopeCount(dimension)));
  }

  FaceIterator MPIMesh::getBoundary() const
  {
    std::vector<Index> indices;
    const auto& shard = getShard();
    const size_t d = getDimension() - 1;
    const size_t count = shard.getFaceCount();
    for (Index i = 0; i < count; ++i)
    {
      if (shard.isOwned(d, i) && shard.isBoundary(i))
        indices.push_back(i);
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  FaceIterator MPIMesh::getInterface() const
  {
    std::vector<Index> indices;
    const auto& shard = getShard();
    const size_t d = getDimension() - 1;
    const size_t count = shard.getFaceCount();
    for (Index i = 0; i < count; ++i)
    {
      if (shard.isOwned(d, i) && shard.isInterface(i))
        indices.push_back(i);
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  bool MPIMesh::isBoundary(Index faceIdx) const
  {
    const auto& shard = getShard();
    const size_t d = getDimension() - 1;
    assert(faceIdx < shard.getPolytopeCount(d));
    return shard.isOwned(d, faceIdx) && shard.isBoundary(faceIdx);
  }

  bool MPIMesh::isInterface(Index faceIdx) const
  {
    const auto& shard = getShard();
    const size_t d = getDimension() - 1;
    assert(faceIdx < shard.getPolytopeCount(d));
    return shard.isOwned(d, faceIdx) && shard.isInterface(faceIdx);
  }

  Real MPIMesh::getVolume() const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = this->getShard();
    Real local = 0;
    for (auto it = shard.getPolytope(3); it; ++it)
    {
      if (shard.isOwned(3, it->getIndex()))
        local += it->getMeasure();
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getVolume(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = this->getShard();
    Real local = 0;
    for (auto it = shard.getPolytope(3); it; ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (shard.isOwned(3, it->getIndex()))
          local += it->getMeasure();
      }
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getVolume(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm  = m_context.getCommunicator();
    const auto& shard = this->getShard();

    Real local = 0;
    for (auto it = shard.getPolytope(3); it; ++it)
    {
      const Index i = it->getIndex();
      if (!shard.isOwned(3, i))
        continue;

      const Optional<Attribute> oa = shard.getAttribute(3, i); // Optional<Attribute>
      if (oa && attrs.contains(*oa))
        local += it->getMeasure();
    }

    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getPerimeter() const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = 0;
    for (auto it = this->getBoundary(); it; ++it)
      local += it->getMeasure();
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getPerimeter(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = 0;
    for (auto it = this->getBoundary(); it; ++it)
    {
      if (it->getAttribute() == attr)
        local += it->getMeasure();
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getPerimeter(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = 0;
    for (auto it = this->getBoundary(); it; ++it)
    {
      const Optional<Attribute> a = it->getAttribute();
      if (a && attrs.contains(*a))
        local += it->getMeasure();
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getArea() const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = this->getShard();
    Real local = 0;
    for (auto it = shard.getPolytope(2); it; ++it)
    {
      if (shard.isOwned(2, it->getIndex()))
        local += it->getMeasure();
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getArea(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = this->getShard();
    Real local = 0;
    for (auto it = shard.getPolytope(2); it; ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (shard.isOwned(2, it->getIndex()))
          local += it->getMeasure();
      }
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getArea(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm  = m_context.getCommunicator();
    const auto& shard = this->getShard();

    Real local = 0;
    for (auto it = shard.getPolytope(2); it; ++it)
    {
      const Index i = it->getIndex();
      if (!shard.isOwned(2, i))
        continue;

      const Optional<Attribute> a = it->getAttribute(); // Optional<Attribute>
      if (a && attrs.contains(*a))
        local += it->getMeasure();
    }

    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getMeasure(size_t d) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = this->getShard();
    Real local = 0;
    for (auto it = shard.getPolytope(d); it; ++it)
    {
      if (shard.isOwned(d, it->getIndex()))
        local += it->getMeasure();
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getMeasure(size_t d, Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = this->getShard();
    Real local = 0;
    for (auto it = shard.getPolytope(d); it; ++it)
    {
      if (it->getAttribute() == attr)
      {
        if (shard.isOwned(d, it->getIndex()))
          local += it->getMeasure();
      }
    }
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getMeasure(size_t d, const FlatSet<Attribute>& attrs) const
  {
    const auto& comm  = m_context.getCommunicator();
    const auto& shard = this->getShard();

    Real local = 0;
    for (auto it = shard.getPolytope(d); it; ++it)
    {
      const Index i = it->getIndex();
      if (!shard.isOwned(d, i))
        continue;

      const Optional<Attribute> a = it->getAttribute(); // Optional<Attribute>
      if (a && attrs.contains(*a))
        local += it->getMeasure();
    }

    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  const PolytopeTransformation& MPIMesh::getPolytopeTransformation(
      size_t dimension, Index localIdx) const
  {
    const auto& shard = this->getShard();
    assert(localIdx < shard.getPolytopeCount(dimension));
    return shard.getPolytopeTransformation(dimension, localIdx);
  }

  Polytope::Type MPIMesh::getGeometry(size_t dimension, Index localIdx) const
  {
    const auto& shard = this->getShard();
    assert(localIdx < shard.getPolytopeCount(dimension));
    return shard.getGeometry(dimension, localIdx);
  }

  Optional<Attribute> MPIMesh::getAttribute(size_t dimension, Index localIdx) const
  {
    const auto& shard = this->getShard();
    assert(localIdx < shard.getPolytopeCount(dimension));
    return shard.getAttribute(dimension, localIdx);
  }

  MPIMesh& MPIMesh::setAttribute(
      const std::pair<size_t, Index>& p,
      const Optional<Attribute>& attr)
  {
    const auto& [d, localIdx] = p;
    auto& shard = this->getShard();
    assert(localIdx < shard.getPolytopeCount(d));
    shard.setAttribute({ d, localIdx }, attr);
    return *this;
  }

  MPIMesh& MPIMesh::setVertexCoordinates(Index localIdx, Real s, size_t i)
  {
    auto& shard = this->getShard();
    assert(localIdx < shard.getVertexCount());
    shard.setVertexCoordinates(localIdx, s, i);
    return *this;
  }

  MPIMesh& MPIMesh::setVertexCoordinates(Index localIdx, const Math::SpatialPoint& coords)
  {
    auto& shard = this->getShard();
    assert(localIdx < shard.getVertexCount());
    shard.setVertexCoordinates(localIdx, coords);
    return *this;
  }

  MPIMesh& MPIMesh::setPolytopeTransformation(
      const std::pair<size_t, Index> p, PolytopeTransformation* trans)
  {
    const auto& [d, localIdx] = p;
    auto& shard = this->getShard();
    assert(localIdx < shard.getPolytopeCount(d));
    shard.setPolytopeTransformation({ d, localIdx }, trans);
    return *this;
  }

  MPIMesh& MPIMesh::load(const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    auto& shard = this->getShard();
    shard.load(filename, fmt);
    return *this;
  }

  void MPIMesh::save(const boost::filesystem::path& filename, IO::FileFormat fmt) const
  {
    const auto& shard = getShard();
    shard.save(filename, fmt);
  }

  Math::SpatialPoint MPIMesh::getVertexCoordinates(Index localIdx) const
  {
    const auto& shard = this->getShard();
    assert(localIdx < shard.getVertexCount());
    return shard.getVertexCoordinates(localIdx);
  }

  Connectivity<Context::Local>& MPIMesh::getConnectivity()
  {
    return m_shard.getConnectivity();
  }

  const Connectivity<Context::Local>& MPIMesh::getConnectivity() const
  {
    return m_shard.getConnectivity();
  }

  Mesh<Context::MPI>& MPIMesh::reconcile(size_t d)
  {
    auto& shard = this->getShard();
    auto& conn  = shard.getConnectivity();

    const auto& ctx  = m_context;
    const auto& comm = ctx.getCommunicator();

    const int rank = comm.rank();
    const size_t D = shard.getDimension();

    comm.barrier();

    assert(d <= D);

    // Vertices and cells are assumed to have already been reconciled by sharding.
    if (d == 0 || d == D)
      return *this;

    const auto& D2d = conn.getIncidence(D, d);
    const auto& d20 = conn.getIncidence(d, 0);

    const size_t nd = shard.getPolytopeCount(d);
    const size_t nc = shard.getPolytopeCount(D);

    if (nd == 0)
      return *this;

    // --------------------------------------------------------------------------
    // Build cell communication stencil from the existing shard metadata.
    // --------------------------------------------------------------------------
    UnorderedSet<int> neighbors;

    {
      const auto& cellHalo  = shard.getHalo(D);   // owned cell -> remote ranks
      const auto& cellOwner = shard.getOwner(D);  // ghost cell -> owner rank

      for (const auto& [cell, rs] : cellHalo)
      {
        (void) cell;
        for (const auto r : rs)
          if (static_cast<int>(r) != rank)
            neighbors.insert(static_cast<int>(r));
      }

      for (const auto& [cell, r] : cellOwner)
      {
        (void) cell;
        if (static_cast<int>(r) != rank)
          neighbors.insert(static_cast<int>(r));
      }
    }

    // --------------------------------------------------------------------------
    // Per-entity local classification from local cell incidences.
    //
    // ownedIncident[i] : entity i touches at least one owned cell
    // ghostIncident[i] : entity i touches at least one ghost cell
    // minGhostOwner[i] : minimum owner rank among incident ghost cells
    // --------------------------------------------------------------------------
    std::vector<uint8_t> ownedIncident(nd, 0);
    std::vector<uint8_t> ghostIncident(nd, 0);
    std::vector<int> minGhostOwner(nd, std::numeric_limits<int>::max());

    for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
    {
      const bool owned = shard.isOwned(D, cell);
      const bool ghost = shard.isGhost(D, cell);

      int gOwner = std::numeric_limits<int>::max();
      if (ghost)
      {
        const auto it = shard.getOwner(D).find(cell);
        assert(it != shard.getOwner(D).end());
        gOwner = static_cast<int>(it->second);
      }

      for (const Index e : D2d[cell])
      {
        if (owned)
          ownedIncident[e] = 1;

        if (ghost)
        {
          ghostIncident[e] = 1;
          minGhostOwner[e] = std::min(minGhostOwner[e], gOwner);
        }
      }
    }

    // --------------------------------------------------------------------------
    // Canonical key for each local entity of dimension d.
    //
    // The key is:
    //   [ geometry_code, sorted(distributed_vertex_ids...) ]
    //
    // This assumes vertex distributed ids are already valid.
    // --------------------------------------------------------------------------
    std::vector<std::vector<Index>> keys(nd);
    UnorderedMap<std::vector<Index>, Index> key2local;

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      const auto g = shard.getGeometry(d, i);
      const auto& lv = d20[i];

      std::vector<Index> key;
      key.reserve(lv.size() + 1);
      key.push_back(static_cast<Index>(g));

      for (const Index localVertex : lv)
      {
        const Index distributedVertex = shard.getPolytopeMap(0).left.at(localVertex);
        key.push_back(distributedVertex);
      }

      std::sort(key.begin() + 1, key.end());

      keys[i] = key;
      key2local.emplace(keys[i], i);
    }

    // --------------------------------------------------------------------------
    // Candidate entities for distributed exchange:
    // anything touching at least one ghost cell.
    //
    // Internal entities (owned-only, no ghost incidence) are purely local.
    // --------------------------------------------------------------------------
    std::vector<std::vector<Index>> localCandidateKeys;
    localCandidateKeys.reserve(nd);

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      if (ghostIncident[i])
        localCandidateKeys.push_back(keys[i]);
    }

    // Sharer sets: ranks known to contain the same distributed entity.
    std::vector<UnorderedSet<int>> sharers(nd);

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      if (ownedIncident[i])
        sharers[i].insert(rank);
    }

    // --------------------------------------------------------------------------
    // Phase 1: exchange candidate keys with neighbors and discover sharers.
    // --------------------------------------------------------------------------
    const int tagKeys = 1000 + static_cast<int>(d);
    {
      std::vector<boost::mpi::request> sends;
      sends.reserve(neighbors.size());

      for (const int nbr : neighbors)
        sends.push_back(comm.isend(nbr, tagKeys, localCandidateKeys));

      for (const int nbr : neighbors)
      {
        std::vector<std::vector<Index>> remoteKeys;
        comm.recv(nbr, tagKeys, remoteKeys);

        for (const auto& key : remoteKeys)
        {
          const auto it = key2local.find(key);
          if (it != key2local.end())
            sharers[it->second].insert(nbr);
        }
      }

      boost::mpi::wait_all(sends.begin(), sends.end());
    }

    // --------------------------------------------------------------------------
    // Decide owner rank for each local entity.
    //
    // Rules:
    // - internal owned-only entity: owner = self
    // - shared entity: owner = min(sharers U ghost-cell-owners)
    // - ghost-only entity: owner = min incident ghost-cell-owner, refined by matches
    // --------------------------------------------------------------------------
    std::vector<int> ownerRank(nd, std::numeric_limits<int>::max());

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      int o = std::numeric_limits<int>::max();

      if (ownedIncident[i])
        o = std::min(o, rank);

      if (minGhostOwner[i] != std::numeric_limits<int>::max())
        o = std::min(o, minGhostOwner[i]);

      for (const int r : sharers[i])
        o = std::min(o, r);

      assert(o != std::numeric_limits<int>::max());
      ownerRank[i] = o;
    }

    // --------------------------------------------------------------------------
    // Assign distributed ids to locally owned entities.
    //
    // Ownership now means: ownerRank[i] == rank.
    // --------------------------------------------------------------------------
    std::vector<Index> distId(nd, std::numeric_limits<Index>::max());
    std::vector<Index> locallyOwned;
    locallyOwned.reserve(nd);

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      if (ownerRank[i] == rank)
        locallyOwned.push_back(i);
    }

    const size_t localOwnedCount = locallyOwned.size();
    std::vector<size_t> ownedCounts;
    boost::mpi::all_gather(comm, localOwnedCount, ownedCounts);

    size_t offset = 0;
    for (int r = 0; r < rank; ++r)
      offset += ownedCounts[r];

    for (size_t k = 0; k < locallyOwned.size(); ++k)
      distId[locallyOwned[k]] = static_cast<Index>(offset + k);

    // --------------------------------------------------------------------------
    // Phase 2: owners send back distributed ids to neighbors that also contain
    // the entity.
    // --------------------------------------------------------------------------
    const int tagIds = 2000 + static_cast<int>(d);
    {
      UnorderedMap<int, std::vector<std::pair<std::vector<Index>, Index>>> sendIds;

      for (Index i = 0; i < static_cast<Index>(nd); ++i)
      {
        if (ownerRank[i] != rank)
          continue;

        for (const int r : sharers[i])
        {
          if (r == rank)
            continue;
          sendIds[r].push_back({ keys[i], distId[i] });
        }
      }

      std::vector<boost::mpi::request> sends;
      sends.reserve(sendIds.size());

      for (const auto& [nbr, payload] : sendIds)
        sends.push_back(comm.isend(nbr, tagIds, payload));

      for (const int nbr : neighbors)
      {
        std::vector<std::pair<std::vector<Index>, Index>> remoteIds;
        comm.recv(nbr, tagIds, remoteIds);

        for (const auto& [key, gid] : remoteIds)
        {
          const auto it = key2local.find(key);
          if (it != key2local.end())
          {
            const Index i = it->second;
            if (ownerRank[i] != rank)
            {
              distId[i] = gid;
              ownerRank[i] = nbr;
            }
          }
        }
      }

      boost::mpi::wait_all(sends.begin(), sends.end());
    }

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
      assert(distId[i] != std::numeric_limits<Index>::max());

    // --------------------------------------------------------------------------
    // Rewrite shard metadata for dimension d.
    // --------------------------------------------------------------------------
    {
      // distributed-index map
      auto& map = shard.getPolytopeMap(d);
      map.left.resize(nd);
      map.right.clear();
      map.right.reserve(nd);

      for (Index i = 0; i < static_cast<Index>(nd); ++i)
      {
        map.left[i] = distId[i];
        const auto [it, inserted] = map.right.insert({ distId[i], i });
        (void) it;
        assert(inserted);
      }

      // flags
      auto& flags = shard.getFlags(d);
      flags.assign(nd, Shard::Flags::None);

      // owner / halo
      auto& owner = shard.getOwner(d);
      auto& halo  = shard.getHalo(d);
      owner.clear();
      halo.clear();

      for (Index i = 0; i < static_cast<Index>(nd); ++i)
      {
        if (ownerRank[i] == rank)
        {
          flags[i] = Shard::Flags::Owned;

          IndexSet rs;
          for (const int r : sharers[i])
          {
            if (r != rank)
              rs.insert(static_cast<Index>(r));
          }

          if (!rs.empty())
            halo.emplace(i, std::move(rs));
        }
        else
        {
          flags[i] = Shard::Flags::Ghost;
          owner.emplace(i, static_cast<Index>(ownerRank[i]));
        }
      }
    }

    return *this;
  }
}

