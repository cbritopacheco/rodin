#include <algorithm>

#include <boost/dynamic_bitset.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Alert/NamespacedException.h"
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
    const auto& shard = this->getShard();
    return CellIterator(*this, BoundedIndexGenerator(localIdx, shard.getCellCount()));
  }

  FaceIterator MPIMesh::getFace(Index localIdx) const
  {
    const auto& shard = this->getShard();
    return FaceIterator(*this, BoundedIndexGenerator(localIdx, shard.getFaceCount()));
  }

  VertexIterator MPIMesh::getVertex(Index localIdx) const
  {
    const auto& shard = this->getShard();
    return VertexIterator(*this, BoundedIndexGenerator(localIdx, shard.getVertexCount()));
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index localIdx) const
  {
    const auto& shard = this->getShard();
    return PolytopeIterator(
        dimension, *this, BoundedIndexGenerator(localIdx, shard.getPolytopeCount(dimension)));
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

  const PolytopeQuadrature&
  MPIMesh::getQuadrature(
      size_t dimension, Index localIdx, const QF::QuadratureFormulaBase& qf) const
  {
    const auto& shard = this->getShard();
    assert(localIdx < shard.getPolytopeCount(dimension));
    return shard.getQuadrature(dimension, localIdx, qf);
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

  Mesh<Context::MPI>& MPIMesh::reconcile(size_t d, const ReconcileOptions& options)
  {
    auto& shard = this->getShard();
    auto& conn  = shard.getConnectivity();

    const auto& comm = m_context.getCommunicator();
    const int rank   = comm.rank();

    const size_t D = shard.getDimension();
    assert(d <= D);

    // Vertices and top cells are assumed already reconciled.
    if (d == 0 || d == D)
      return *this;

    // Need:
    // - D -> d to visit subentities through incident top cells
    // - d -> 0 to build canonical keys from distributed vertex ids
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(shard, D, d);
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(shard, d, 0);

    const auto& D2d = conn.getIncidence(D, d);
    const auto& d20 = conn.getIncidence(d, 0);

    const size_t nd = shard.getPolytopeCount(d);
    const size_t nc = shard.getPolytopeCount(D);

    if (nd == 0)
      return *this;

    using Key = Polytope::Key;
    using KeyMap = UnorderedMap<
        Key, Index,
        Polytope::Key::SymmetricHash,
        Polytope::Key::SymmetricEquality>;

    static constexpr int kInfOwner = std::numeric_limits<int>::max();
    static constexpr Index kInvalidId = std::numeric_limits<Index>::max();

    // --------------------------------------------------------------------------
    // Neighbor stencil from top-cell ownership metadata.
    // --------------------------------------------------------------------------
    std::vector<int> neighbors;
    UnorderedMap<int, size_t> neighborPos;
    {
      UnorderedSet<int> nbrs;

      const auto& cellHalo  = shard.getHalo(D);   // owned top cell -> remote ranks that also store it
      const auto& cellOwner = shard.getOwner(D);  // ghost top cell -> owner rank

      for (const auto& [cell, rs] : cellHalo)
      {
        (void) cell;
        for (const Index r : rs)
        {
          const int rr = static_cast<int>(r);
          if (rr != rank)
            nbrs.insert(rr);
        }
      }

      for (const auto& [cell, r] : cellOwner)
      {
        (void) cell;
        const int rr = static_cast<int>(r);
        if (rr != rank)
          nbrs.insert(rr);
      }

      neighbors.assign(nbrs.begin(), nbrs.end());
      std::sort(neighbors.begin(), neighbors.end());

      neighborPos.reserve(neighbors.size());
      for (size_t k = 0; k < neighbors.size(); ++k)
        neighborPos.emplace(neighbors[k], k);
    }

    // --------------------------------------------------------------------------
    // Canonical key of each local d-entity from distributed vertex ids.
    // --------------------------------------------------------------------------
    std::vector<Key> keys(nd);
    KeyMap key2local;
    key2local.reserve(nd);

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      const auto& lv = d20[i];
      Key key(static_cast<std::uint8_t>(lv.size()));

      for (std::uint8_t j = 0; j < static_cast<std::uint8_t>(lv.size()); ++j)
      {
        const Index localVertex = lv[j];
        key[j] = shard.getPolytopeMap(0).left.at(localVertex);
      }

      keys[i] = key;
      const auto [it, inserted] = key2local.emplace(keys[i], i);
      (void) it;
      assert(inserted);
    }

    // --------------------------------------------------------------------------
    // Per-entity local classification and per-neighbor send lists.
    //
    // inLocalPartition[i]:
    //   true iff entity i is incident to at least one locally owned top cell.
    //
    // entitiesToSend[k]:
    //   deduplicated list of local entity indices to export to neighbors[k].
    //   Built from the cell-entity-neighbor stencil so that each entity-
    //   neighbor pair appears exactly once.
    // --------------------------------------------------------------------------
    std::vector<uint8_t> inLocalPartition(nd, 0);
    std::vector<std::vector<Index>> entitiesToSend(neighbors.size());

    {
      const auto& cellHalo  = shard.getHalo(D);
      const auto& cellOwner = shard.getOwner(D);

      for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
      {
        const bool owned = shard.isOwned(D, cell);
        const bool ghost = shard.isGhost(D, cell);

        if (owned)
        {
          for (const Index e : D2d[cell])
            inLocalPartition[e] = 1;

          auto hit = cellHalo.find(cell);
          if (hit != cellHalo.end())
          {
            for (const Index r : hit->second)
            {
              const int rr = static_cast<int>(r);
              if (rr == rank)
                continue;
              const auto nit = neighborPos.find(rr);
              assert(nit != neighborPos.end());
              for (const Index e : D2d[cell])
                entitiesToSend[nit->second].push_back(e);
            }
          }
        }

        if (ghost)
        {
          auto oit = cellOwner.find(cell);
          assert(oit != cellOwner.end());
          const int rr = static_cast<int>(oit->second);
          if (rr != rank)
          {
            const auto nit = neighborPos.find(rr);
            assert(nit != neighborPos.end());
            for (const Index e : D2d[cell])
              entitiesToSend[nit->second].push_back(e);
          }
        }
      }

      for (auto& v : entitiesToSend)
      {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
      }
    }

    // --------------------------------------------------------------------------
    // Key exchange with candidate neighbors.
    //
    // Each rank sends (key, inLocalPartition flag, senderLocalIdx) for entities
    // in its send list.  Recipients match received keys against their local key
    // map and build compact per-neighbor structures for subsequent rounds.
    //
    // convergeSend[k]:
    //   list of local entity indices confirmed shared with neighbors[k].
    //   Replaces per-entity matchedHolders sets with K per-neighbor vectors.
    //
    // remoteToLocal[k]:
    //   map from neighbor k's sender local index to our local entity index.
    //   Enables index-based (not key-based) addressing in convergence rounds,
    //   reducing per-entry message size from ~69 bytes to ~12 bytes.
    //
    // ownerRank[i]:
    //   seeded from local inLocalPartition and remote flags received here.
    //   For most entities this already yields the correct owner; iterative
    //   convergence below refines the remaining multi-hop cases.
    // --------------------------------------------------------------------------
    // KeyMsg: (canonical key, (inLocalPartition flag, sender's local index))
    using KeyMsg = std::pair<Key, std::pair<uint8_t, Index>>;

    std::vector<std::vector<Index>> convergeSend(neighbors.size());
    std::vector<UnorderedMap<Index, Index>> remoteToLocal(neighbors.size());

    std::vector<int> ownerRank(nd, kInfOwner);
    for (Index i = 0; i < static_cast<Index>(nd); ++i)
      ownerRank[i] = inLocalPartition[i] ? rank : kInfOwner;

    {
      std::vector<std::vector<KeyMsg>> sendbuf(neighbors.size());
      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        sendbuf[k].reserve(entitiesToSend[k].size());
        for (const Index i : entitiesToSend[k])
          sendbuf[k].push_back({ keys[i], { inLocalPartition[i], i } });
      }

      const int tagKeys = 1000 + static_cast<int>(d);

      std::vector<boost::mpi::request> reqs;
      reqs.reserve(neighbors.size());

      for (size_t k = 0; k < neighbors.size(); ++k)
        reqs.push_back(comm.isend(neighbors[k], tagKeys, sendbuf[k]));

      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        std::vector<KeyMsg> recvbuf;
        comm.recv(neighbors[k], tagKeys, recvbuf);

        const int nbr = neighbors[k];
        for (const auto& [key, flagAndIdx] : recvbuf)
        {
          const auto it = key2local.find(key);
          if (it == key2local.end())
            continue;

          const Index i = it->second;
          convergeSend[k].push_back(i);
          remoteToLocal[k].emplace(flagAndIdx.second, i);

          if (flagAndIdx.first && nbr < ownerRank[i])
            ownerRank[i] = nbr;
        }
      }

      boost::mpi::wait_all(reqs.begin(), reqs.end());
    }

    // --------------------------------------------------------------------------
    // Iterative owner convergence.
    //
    // Repeatedly exchange the current best owner candidate with matched holders
    // and take the global minimum.  This propagates the smallest partition-
    // holder rank across the full holder graph, even when it spans multiple
    // neighbor hops (e.g. edges in 3D shared by a ring of tetrahedra whose
    // partitions are not face-adjacent).
    //
    // Messages use (senderLocalIdx, owner) pairs instead of (Key, owner),
    // reducing per-entry size from ~69 bytes to ~12 bytes.  The receiver
    // resolves senderLocalIdx via the remoteToLocal map built during key
    // exchange.
    //
    // Because ownerRank was already seeded from both local and received flags
    // above, most entities already carry the correct owner.  The loop typically
    // converges in one to two additional rounds.
    // --------------------------------------------------------------------------
    using OwnerMsg = std::pair<Index, int>;

    const size_t checkPeriod = std::max<size_t>(1, options.globalCheckPeriod);
    size_t ownerRounds = 0;
    bool ownerDirtySinceCheck = false;
    while (true)
    {
      if (ownerRounds >= options.maxOwnerRounds)
      {
        if (options.strictRoundCap)
        {
          throw std::runtime_error(
              "MPIMesh::reconcile owner convergence exceeded maxOwnerRounds");
        }
      }
      ++ownerRounds;

      std::vector<std::vector<OwnerMsg>> sendbuf(neighbors.size());

      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        sendbuf[k].reserve(convergeSend[k].size());
        for (const Index i : convergeSend[k])
          sendbuf[k].push_back({ i, ownerRank[i] });
      }

      const int tagOwner = 1500 + static_cast<int>(d);

      std::vector<boost::mpi::request> reqs;
      reqs.reserve(neighbors.size());

      for (size_t k = 0; k < neighbors.size(); ++k)
        reqs.push_back(comm.isend(neighbors[k], tagOwner, sendbuf[k]));

      bool changed = false;

      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        std::vector<OwnerMsg> recvbuf;
        comm.recv(neighbors[k], tagOwner, recvbuf);

        for (const auto& [remoteIdx, remoteOwner] : recvbuf)
        {
          const auto it = remoteToLocal[k].find(remoteIdx);
          if (it == remoteToLocal[k].end())
            continue;

          const Index i = it->second;
          if (remoteOwner < ownerRank[i])
          {
            ownerRank[i] = remoteOwner;
            changed = true;
          }
        }
      }

      boost::mpi::wait_all(reqs.begin(), reqs.end());

      ownerDirtySinceCheck = ownerDirtySinceCheck || changed;

      const bool forceCheck = !changed;
      const bool periodicCheck = (ownerRounds % checkPeriod) == 0;
      if (forceCheck || periodicCheck)
      {
        bool globallyChanged = false;
        boost::mpi::all_reduce(comm, ownerDirtySinceCheck, globallyChanged, std::logical_or<bool>());
        ownerDirtySinceCheck = false;

        if (!globallyChanged)
          break;
      }
    }

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
      assert(ownerRank[i] != kInfOwner);

    // --------------------------------------------------------------------------
    // Assign distributed ids on owners only.
    // --------------------------------------------------------------------------
    std::vector<Index> distId(nd, kInvalidId);
    std::vector<Index> locallyOwned;
    locallyOwned.reserve(nd);

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      if (ownerRank[i] == rank)
        locallyOwned.push_back(i);
    }

    const size_t localOwnedCount = locallyOwned.size();

    size_t inclusive = 0;
    boost::mpi::scan(comm, localOwnedCount, inclusive, std::plus<size_t>());

    const size_t offset = inclusive - localOwnedCount;

    for (size_t k = 0; k < locallyOwned.size(); ++k)
      distId[locallyOwned[k]] = static_cast<Index>(offset + k);

    // --------------------------------------------------------------------------
    // Iterative gid convergence.
    //
    // Owners start with a valid gid. Repeatedly exchange gids with matched
    // holders until every local copy has received the gid.  Uses index-based
    // messages (senderLocalIdx, gid) for the same bandwidth savings as above.
    // --------------------------------------------------------------------------
    using GidMsg = std::pair<Index, Index>;

    size_t gidRounds = 0;
    bool gidDirtySinceCheck = false;
    while (true)
    {
      if (gidRounds >= options.maxGidRounds)
      {
        if (options.strictRoundCap)
        {
          throw std::runtime_error(
              "MPIMesh::reconcile gid convergence exceeded maxGidRounds");
        }
      }
      ++gidRounds;

      std::vector<std::vector<GidMsg>> sendbuf(neighbors.size());

      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        for (const Index i : convergeSend[k])
        {
          if (distId[i] != kInvalidId)
            sendbuf[k].push_back({ i, distId[i] });
        }
      }

      const int tagIds = 2000 + static_cast<int>(d);

      std::vector<boost::mpi::request> reqs;
      reqs.reserve(neighbors.size());

      for (size_t k = 0; k < neighbors.size(); ++k)
        reqs.push_back(comm.isend(neighbors[k], tagIds, sendbuf[k]));

      bool changed = false;

      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        std::vector<GidMsg> recvbuf;
        comm.recv(neighbors[k], tagIds, recvbuf);

        for (const auto& [remoteIdx, gid] : recvbuf)
        {
          const auto it = remoteToLocal[k].find(remoteIdx);
          if (it == remoteToLocal[k].end())
            continue;

          const Index i = it->second;
          if (distId[i] == kInvalidId)
          {
            distId[i] = gid;
            changed = true;
          }
          else
          {
            assert(distId[i] == gid);
          }
        }
      }

      boost::mpi::wait_all(reqs.begin(), reqs.end());

      gidDirtySinceCheck = gidDirtySinceCheck || changed;

      const bool forceCheck = !changed;
      const bool periodicCheck = (gidRounds % checkPeriod) == 0;
      if (forceCheck || periodicCheck)
      {
        bool globallyChanged = false;
        boost::mpi::all_reduce(comm, gidDirtySinceCheck, globallyChanged, std::logical_or<bool>());
        gidDirtySinceCheck = false;

        if (!globallyChanged)
          break;
      }
    }

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
      assert(distId[i] != kInvalidId);

    // --------------------------------------------------------------------------
    // Rewrite shard metadata for dimension d.
    // --------------------------------------------------------------------------
    {
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

      auto& state = shard.getState(d);
      state.assign(nd, Shard::State::Ghost);

      auto& owner = shard.getOwner(d);
      auto& halo  = shard.getHalo(d);
      owner.clear();
      halo.clear();

      for (Index i = 0; i < static_cast<Index>(nd); ++i)
      {
        if (ownerRank[i] == rank)
        {
          assert(inLocalPartition[i]);
          state[i] = Shard::State::Owned;
        }
        else
        {
          state[i] = inLocalPartition[i]
            ? Shard::State::Shared
            : Shard::State::Ghost;

          owner.emplace(i, static_cast<Index>(ownerRank[i]));
        }
      }

      // Build halo only for owned entities from convergeSend.
      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        const Index nbr = static_cast<Index>(neighbors[k]);
        for (const Index i : convergeSend[k])
        {
          if (ownerRank[i] == rank)
            halo[i].insert(nbr);
        }
      }
    }

    return *this;
  }

  MPIMesh MPIMesh::UniformGrid(
      const Context::MPI& context, Polytope::Type g, const Array<size_t>& shape)
  {
    const auto& comm = context.getCommunicator();
    const int rank = comm.rank();
    const int size = comm.size();

    struct Interval
    {
      size_t begin = 0;
      size_t end   = 0; // half-open [begin, end)

      size_t size() const
      {
        return end >= begin ? end - begin : 0;
      }

      bool contains(size_t i) const
      {
        return begin <= i && i < end;
      }

      bool empty() const
      {
        return begin == end;
      }
    };

    const auto makeIndexArray = [](auto... xs)
    {
      IndexArray res(sizeof...(xs));
      const Index data[] = { static_cast<Index>(xs)... };
      for (size_t i = 0; i < sizeof...(xs); ++i)
        res(static_cast<Eigen::Index>(i)) = data[i];
      return res;
    };

    const auto finish = [&](Shard&& shard) -> MPIMesh
    {
      return MPIMesh::Builder(context)
        .initialize(std::move(shard))
        .finalize();
    };

    const auto splitRange = [](size_t n, int r, int p) -> Interval
    {
      const size_t q   = n / static_cast<size_t>(p);
      const size_t rem = n % static_cast<size_t>(p);
      const size_t b   = static_cast<size_t>(r) * q + std::min(static_cast<size_t>(r), rem);
      const size_t c   = q + (static_cast<size_t>(r) < rem ? 1 : 0);
      return Interval{ b, b + c };
    };

    const auto expandByOne = [](const Interval& I, size_t n) -> Interval
    {
      if (I.empty())
        return I;
      return Interval{
        I.begin > 0 ? I.begin - 1 : 0,
        I.end   < n ? I.end + 1   : n
      };
    };

    const auto factorize = [](int n)
    {
      std::vector<int> f;
      for (int p = 2; p * p <= n; ++p)
      {
        while (n % p == 0)
        {
          f.push_back(p);
          n /= p;
        }
      }
      if (n > 1)
        f.push_back(n);
      std::sort(f.rbegin(), f.rend());
      return f;
    };

    const auto chooseProcShape = [&](size_t dim, const std::vector<size_t>& cells)
    {
      std::vector<int> ps(dim, 1);
      auto factors = factorize(size);

      for (int f : factors)
      {
        int best = 0;
        double bestScore = -1.0;

        for (size_t d = 0; d < dim; ++d)
        {
          const double score = static_cast<double>(cells[d]) / static_cast<double>(ps[d]);
          const bool feasible = (static_cast<size_t>(ps[d] * f) <= std::max<size_t>(cells[d], 1));
          const double adjusted = feasible ? score : 0.5 * score;

          if (adjusted > bestScore)
          {
            bestScore = adjusted;
            best = static_cast<int>(d);
          }
        }

        ps[best] *= f;
      }

      return ps;
    };

    const auto rankToProcCoord = [](int r, const std::vector<int>& ps)
    {
      std::vector<int> c(ps.size(), 0);
      for (size_t d = 0; d < ps.size(); ++d)
      {
        c[d] = r % ps[d];
        r /= ps[d];
      }
      return c;
    };

    const auto procCoordToRank = [](const std::vector<int>& c, const std::vector<int>& ps)
    {
      int r = 0;
      int stride = 1;
      for (size_t d = 0; d < ps.size(); ++d)
      {
        r += c[d] * stride;
        stride *= ps[d];
      }
      return r;
    };

    const auto uniqueSort = [](std::vector<int>& v)
    {
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
    };

    const auto sp0 = []()
    {
      return Math::SpatialPoint();
    };

    const auto sp1 = [](Real x0)
    {
      Math::SpatialPoint x(1);
      x(0) = x0;
      return x;
    };

    const auto sp2 = [](Real x0, Real x1)
    {
      Math::SpatialPoint x(2);
      x(0) = x0;
      x(1) = x1;
      return x;
    };

    const auto sp3 = [](Real x0, Real x1, Real x2)
    {
      Math::SpatialPoint x(3);
      x(0) = x0;
      x(1) = x1;
      x(2) = x2;
      return x;
    };

    const auto vid1 = [](size_t i) -> Index
    {
      return static_cast<Index>(i);
    };

    const auto vid2 = [](size_t i, size_t j, size_t nx) -> Index
    {
      return static_cast<Index>(i + j * nx);
    };

    const auto vid3 = [](size_t i, size_t j, size_t k, size_t nx, size_t ny) -> Index
    {
      return static_cast<Index>(i + j * nx + k * nx * ny);
    };

    const auto macroId2 = [](size_t i, size_t j, size_t cx) -> Index
    {
      return static_cast<Index>(i + j * cx);
    };

    const auto macroId3 = [](size_t i, size_t j, size_t k, size_t cx, size_t cy) -> Index
    {
      return static_cast<Index>(i + j * cx + k * cx * cy);
    };

    switch (g)
    {
      case Polytope::Type::Point:
      {
        if (shape.size() != 0)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "Expected 0 dimensions for geometry type " << g
            << ", but got " << shape.size() << "."
            << Alert::Raise;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/0, /*sdim=*/0);

        const Shard::State state =
          (rank == 0 ? Shard::State::Owned : Shard::State::Ghost);

        const Index lv = sb.vertex(/*globalIdx=*/0, sp0(), state);

        if (rank == 0)
        {
          for (int r = 1; r < size; ++r)
            sb.halo(0, lv, static_cast<Index>(r));
        }
        else
        {
          sb.setOwner(0, lv, 0);
        }

        return finish(sb.finalize());
      }

      case Polytope::Type::Segment:
      {
        if (shape.size() != 1)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "Expected 1 dimension for geometry type " << g
            << ", but got " << shape.size() << "."
            << Alert::Raise;
        }

        const size_t nx = shape.coeff(0);
        if (nx < 2)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "Segment uniform grid requires at least 2 vertices."
            << Alert::Raise;
        }

        const size_t dim = 1;
        const std::vector<size_t> nCells = { nx - 1 };
        const std::vector<int> procShape = chooseProcShape(dim, nCells);
        const std::vector<int> procCoord = rankToProcCoord(rank, procShape);

        std::vector<std::vector<Interval>> cellSplits(dim);
        for (size_t d = 0; d < dim; ++d)
        {
          cellSplits[d].resize(procShape[d]);
          for (int p = 0; p < procShape[d]; ++p)
            cellSplits[d][p] = splitRange(nCells[d], p, procShape[d]);
        }

        const Interval ownedCells = cellSplits[0][procCoord[0]];
        const Interval ghostCells = expandByOne(ownedCells, nCells[0]);

        const Interval localVerts = Interval{
          ownedCells.begin,
          ownedCells.empty() ? ownedCells.end : std::min(nx, ownedCells.end + 1)
        };

        const Interval ghostVerts = Interval{
          ghostCells.begin,
          ghostCells.empty() ? ghostCells.end : std::min(nx, ghostCells.end + 1)
        };

        const auto ownerOfCell = [&](size_t i) -> int
        {
          for (int px = 0; px < procShape[0]; ++px)
          {
            if (cellSplits[0][px].contains(i))
              return px;
          }
          assert(false);
          return 0;
        };

        const auto ownerOfVertex = [&](size_t i) -> int
        {
          std::vector<int> rs;
          if (i > 0 && i - 1 < nCells[0])
            rs.push_back(ownerOfCell(i - 1));
          if (i < nCells[0])
            rs.push_back(ownerOfCell(i));
          if (rs.empty())
            rs.push_back(0);
          uniqueSort(rs);
          return rs.front();
        };

        const auto holdersOfVertex = [&](size_t i)
        {
          std::vector<int> res;
          for (int p = 0; p < procShape[0]; ++p)
          {
            const Interval Gc = expandByOne(cellSplits[0][p], nCells[0]);
            const Interval Gv{
              Gc.begin,
              Gc.empty() ? Gc.end : std::min(nx, Gc.end + 1)
            };
            if (Gv.contains(i))
              res.push_back(p);
          }
          uniqueSort(res);
          return res;
        };

        const auto stateOfVertex = [&](size_t i) -> Shard::State
        {
          const int owner = ownerOfVertex(i);
          if (localVerts.contains(i))
            return owner == rank ? Shard::State::Owned : Shard::State::Shared;
          return Shard::State::Ghost;
        };

        const auto holdersOfCell = [&](size_t i)
        {
          std::vector<int> res;
          const int owner = ownerOfCell(i);
          const std::vector<int> oc = rankToProcCoord(owner, procShape);

          for (int dx = -1; dx <= 1; ++dx)
          {
            std::vector<int> pc = oc;
            pc[0] += dx;
            if (pc[0] < 0 || pc[0] >= procShape[0])
              continue;

            const Interval G = expandByOne(cellSplits[0][pc[0]], nCells[0]);
            if (G.contains(i))
              res.push_back(procCoordToRank(pc, procShape));
          }

          uniqueSort(res);
          return res;
        };

        Shard::Builder sb;
        sb.initialize(/*dimension=*/1, /*sdim=*/1);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve(ghostVerts.size());

        for (size_t i = ghostVerts.begin; i < ghostVerts.end; ++i)
        {
          const Index gvid = vid1(i);
          const int owner = ownerOfVertex(i);
          const Shard::State state = stateOfVertex(i);

          const Index lv = sb.vertex(gvid, sp1(static_cast<Real>(i)), state);
          gv2lv.emplace(gvid, lv);

          if (state == Shard::State::Owned)
          {
            for (const int r : holdersOfVertex(i))
            {
              if (r != rank)
                sb.halo(0, lv, static_cast<Index>(r));
            }
          }
          else
          {
            sb.setOwner(0, lv, static_cast<Index>(owner));
          }
        }

        for (size_t i = ghostCells.begin; i < ghostCells.end; ++i)
        {
          const int owner = ownerOfCell(i);
          const Shard::State state =
            (owner == rank ? Shard::State::Owned : Shard::State::Ghost);

          const Index gcid = static_cast<Index>(i);
          const IndexArray vs = makeIndexArray(
            gv2lv.at(vid1(i)),
            gv2lv.at(vid1(i + 1))
          );

          const Index lc = sb.polytope(1, gcid, Polytope::Type::Segment, vs, state);

          if (state == Shard::State::Owned)
          {
            for (const int r : holdersOfCell(i))
            {
              if (r != rank)
                sb.halo(1, lc, static_cast<Index>(r));
            }
          }
          else
          {
            sb.setOwner(1, lc, static_cast<Index>(owner));
          }
        }

        return finish(sb.finalize());
      }

      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
      {
        if (shape.size() != 2)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "Expected 2 dimensions for geometry type " << g
            << ", but got " << shape.size() << "."
            << Alert::Raise;
        }

        const size_t nx = shape.coeff(0);
        const size_t ny = shape.coeff(1);
        if (nx < 2 || ny < 2)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "2D uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;

        const size_t dim = 2;
        const std::vector<size_t> nCells = { cx, cy };
        const std::vector<int> procShape = chooseProcShape(dim, nCells);
        const std::vector<int> procCoord = rankToProcCoord(rank, procShape);

        std::vector<std::vector<Interval>> cellSplits(dim);
        for (size_t d = 0; d < dim; ++d)
        {
          cellSplits[d].resize(procShape[d]);
          for (int p = 0; p < procShape[d]; ++p)
            cellSplits[d][p] = splitRange(nCells[d], p, procShape[d]);
        }

        const Interval ownedX = cellSplits[0][procCoord[0]];
        const Interval ownedY = cellSplits[1][procCoord[1]];
        const Interval ghostX = expandByOne(ownedX, cx);
        const Interval ghostY = expandByOne(ownedY, cy);

        const Interval localVertX{
          ownedX.begin,
          ownedX.empty() ? ownedX.end : std::min(nx, ownedX.end + 1)
        };
        const Interval localVertY{
          ownedY.begin,
          ownedY.empty() ? ownedY.end : std::min(ny, ownedY.end + 1)
        };

        const Interval vertX{
          ghostX.begin,
          ghostX.empty() ? ghostX.end : std::min(nx, ghostX.end + 1)
        };
        const Interval vertY{
          ghostY.begin,
          ghostY.empty() ? ghostY.end : std::min(ny, ghostY.end + 1)
        };

        const auto ownerOfCell = [&](size_t i, size_t j) -> int
        {
          int px = -1, py = -1;
          for (int p = 0; p < procShape[0]; ++p)
            if (cellSplits[0][p].contains(i)) { px = p; break; }
          for (int p = 0; p < procShape[1]; ++p)
            if (cellSplits[1][p].contains(j)) { py = p; break; }
          assert(px >= 0 && py >= 0);
          return procCoordToRank({ px, py }, procShape);
        };

        const auto ownerOfVertex = [&](size_t i, size_t j) -> int
        {
          std::vector<int> rs;
          const size_t ix0 = (i > 0 ? i - 1 : i);
          const size_t ix1 = std::min(i, cx ? cx - 1 : 0);
          const size_t jy0 = (j > 0 ? j - 1 : j);
          const size_t jy1 = std::min(j, cy ? cy - 1 : 0);

          for (size_t ii : { ix0, ix1 })
          {
            if (ii >= cx) continue;
            for (size_t jj : { jy0, jy1 })
            {
              if (jj >= cy) continue;
              rs.push_back(ownerOfCell(ii, jj));
            }
          }

          if (rs.empty())
            rs.push_back(0);
          uniqueSort(rs);
          return rs.front();
        };

        const auto holdersOfVertex = [&](size_t i, size_t j)
        {
          std::vector<int> res;
          for (int px = 0; px < procShape[0]; ++px)
          {
            const Interval Gx = expandByOne(cellSplits[0][px], cx);
            const Interval Vx{
              Gx.begin,
              Gx.empty() ? Gx.end : std::min(nx, Gx.end + 1)
            };
            if (!Vx.contains(i))
              continue;

            for (int py = 0; py < procShape[1]; ++py)
            {
              const Interval Gy = expandByOne(cellSplits[1][py], cy);
              const Interval Vy{
                Gy.begin,
                Gy.empty() ? Gy.end : std::min(ny, Gy.end + 1)
              };
              if (Vy.contains(j))
                res.push_back(procCoordToRank({ px, py }, procShape));
            }
          }
          uniqueSort(res);
          return res;
        };

        const auto stateOfVertex = [&](size_t i, size_t j) -> Shard::State
        {
          const int owner = ownerOfVertex(i, j);
          if (localVertX.contains(i) && localVertY.contains(j))
            return owner == rank ? Shard::State::Owned : Shard::State::Shared;
          return Shard::State::Ghost;
        };

        const auto holdersOfCell = [&](size_t i, size_t j)
        {
          std::vector<int> res;
          const int owner = ownerOfCell(i, j);
          const std::vector<int> oc = rankToProcCoord(owner, procShape);

          for (int dx = -1; dx <= 1; ++dx)
          {
            for (int dy = -1; dy <= 1; ++dy)
            {
              std::vector<int> pc = oc;
              pc[0] += dx;
              pc[1] += dy;
              if (pc[0] < 0 || pc[0] >= procShape[0]) continue;
              if (pc[1] < 0 || pc[1] >= procShape[1]) continue;

              const Interval GX = expandByOne(cellSplits[0][pc[0]], cx);
              const Interval GY = expandByOne(cellSplits[1][pc[1]], cy);
              if (GX.contains(i) && GY.contains(j))
                res.push_back(procCoordToRank(pc, procShape));
            }
          }

          uniqueSort(res);
          return res;
        };

        Shard::Builder sb;
        sb.initialize(/*dimension=*/2, /*sdim=*/2);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve(vertX.size() * vertY.size());

        for (size_t j = vertY.begin; j < vertY.end; ++j)
        {
          for (size_t i = vertX.begin; i < vertX.end; ++i)
          {
            const Index gvid = vid2(i, j, nx);
            const int owner = ownerOfVertex(i, j);
            const Shard::State state = stateOfVertex(i, j);

            const Index lv = sb.vertex(
              gvid,
              sp2(static_cast<Real>(i), static_cast<Real>(j)),
              state);
            gv2lv.emplace(gvid, lv);

            if (state == Shard::State::Owned)
            {
              for (const int r : holdersOfVertex(i, j))
              {
                if (r != rank)
                  sb.halo(0, lv, static_cast<Index>(r));
              }
            }
            else
            {
              sb.setOwner(0, lv, static_cast<Index>(owner));
            }
          }
        }

        for (size_t j = ghostY.begin; j < ghostY.end; ++j)
        {
          for (size_t i = ghostX.begin; i < ghostX.end; ++i)
          {
            const Index mid = macroId2(i, j, cx);
            const int owner = ownerOfCell(i, j);
            const Shard::State state =
              (owner == rank ? Shard::State::Owned : Shard::State::Ghost);

            if (g == Polytope::Type::Triangle)
            {
              const Index v00 = gv2lv.at(vid2(i,     j,     nx));
              const Index v10 = gv2lv.at(vid2(i + 1, j,     nx));
              const Index v01 = gv2lv.at(vid2(i,     j + 1, nx));
              const Index v11 = gv2lv.at(vid2(i + 1, j + 1, nx));

              const Index lc0 = sb.polytope(
                2, static_cast<Index>(2 * mid + 0), Polytope::Type::Triangle,
                makeIndexArray(v00, v10, v01), state);

              const Index lc1 = sb.polytope(
                2, static_cast<Index>(2 * mid + 1), Polytope::Type::Triangle,
                makeIndexArray(v10, v11, v01), state);

              if (state == Shard::State::Owned)
              {
                for (const int r : holdersOfCell(i, j))
                {
                  if (r != rank)
                  {
                    sb.halo(2, lc0, static_cast<Index>(r));
                    sb.halo(2, lc1, static_cast<Index>(r));
                  }
                }
              }
              else
              {
                const Index ow = static_cast<Index>(owner);
                sb.setOwner(2, lc0, ow);
                sb.setOwner(2, lc1, ow);
              }
            }
            else
            {
              const IndexArray vs = makeIndexArray(
                gv2lv.at(vid2(i,     j,     nx)),
                gv2lv.at(vid2(i + 1, j,     nx)),
                gv2lv.at(vid2(i + 1, j + 1, nx)),
                gv2lv.at(vid2(i,     j + 1, nx))
              );

              const Index lc = sb.polytope(2, mid, Polytope::Type::Quadrilateral, vs, state);

              if (state == Shard::State::Owned)
              {
                for (const int r : holdersOfCell(i, j))
                {
                  if (r != rank)
                    sb.halo(2, lc, static_cast<Index>(r));
                }
              }
              else
              {
                sb.setOwner(2, lc, static_cast<Index>(owner));
              }
            }
          }
        }

        return finish(sb.finalize());
      }

      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
      {
        if (shape.size() != 3)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "Expected 3 dimensions for geometry type " << g
            << ", but got " << shape.size() << "."
            << Alert::Raise;
        }

        const size_t nx = shape.coeff(0);
        const size_t ny = shape.coeff(1);
        const size_t nz = shape.coeff(2);
        if (nx < 2 || ny < 2 || nz < 2)
        {
          Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
            << "3D uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;
        const size_t cz = nz - 1;

        const size_t dim = 3;
        const std::vector<size_t> nCells = { cx, cy, cz };
        const std::vector<int> procShape = chooseProcShape(dim, nCells);
        const std::vector<int> procCoord = rankToProcCoord(rank, procShape);

        std::vector<std::vector<Interval>> cellSplits(dim);
        for (size_t d = 0; d < dim; ++d)
        {
          cellSplits[d].resize(procShape[d]);
          for (int p = 0; p < procShape[d]; ++p)
            cellSplits[d][p] = splitRange(nCells[d], p, procShape[d]);
        }

        const Interval ownedX = cellSplits[0][procCoord[0]];
        const Interval ownedY = cellSplits[1][procCoord[1]];
        const Interval ownedZ = cellSplits[2][procCoord[2]];

        const Interval ghostX = expandByOne(ownedX, cx);
        const Interval ghostY = expandByOne(ownedY, cy);
        const Interval ghostZ = expandByOne(ownedZ, cz);

        const Interval localVertX{
          ownedX.begin,
          ownedX.empty() ? ownedX.end : std::min(nx, ownedX.end + 1)
        };
        const Interval localVertY{
          ownedY.begin,
          ownedY.empty() ? ownedY.end : std::min(ny, ownedY.end + 1)
        };
        const Interval localVertZ{
          ownedZ.begin,
          ownedZ.empty() ? ownedZ.end : std::min(nz, ownedZ.end + 1)
        };

        const Interval vertX{
          ghostX.begin,
          ghostX.empty() ? ghostX.end : std::min(nx, ghostX.end + 1)
        };
        const Interval vertY{
          ghostY.begin,
          ghostY.empty() ? ghostY.end : std::min(ny, ghostY.end + 1)
        };
        const Interval vertZ{
          ghostZ.begin,
          ghostZ.empty() ? ghostZ.end : std::min(nz, ghostZ.end + 1)
        };

        const auto ownerOfCell = [&](size_t i, size_t j, size_t k) -> int
        {
          int px = -1, py = -1, pz = -1;
          for (int p = 0; p < procShape[0]; ++p)
            if (cellSplits[0][p].contains(i)) { px = p; break; }
          for (int p = 0; p < procShape[1]; ++p)
            if (cellSplits[1][p].contains(j)) { py = p; break; }
          for (int p = 0; p < procShape[2]; ++p)
            if (cellSplits[2][p].contains(k)) { pz = p; break; }
          assert(px >= 0 && py >= 0 && pz >= 0);
          return procCoordToRank({ px, py, pz }, procShape);
        };

        const auto ownerOfVertex = [&](size_t i, size_t j, size_t k) -> int
        {
          std::vector<int> rs;

          std::array<size_t, 2> ii = { i > 0 ? i - 1 : i, std::min(i, cx ? cx - 1 : 0) };
          std::array<size_t, 2> jj = { j > 0 ? j - 1 : j, std::min(j, cy ? cy - 1 : 0) };
          std::array<size_t, 2> kk = { k > 0 ? k - 1 : k, std::min(k, cz ? cz - 1 : 0) };

          for (size_t a : ii)
          {
            if (a >= cx) continue;
            for (size_t b : jj)
            {
              if (b >= cy) continue;
              for (size_t c : kk)
              {
                if (c >= cz) continue;
                rs.push_back(ownerOfCell(a, b, c));
              }
            }
          }

          if (rs.empty())
            rs.push_back(0);
          uniqueSort(rs);
          return rs.front();
        };

        const auto holdersOfVertex = [&](size_t i, size_t j, size_t k)
        {
          std::vector<int> res;
          for (int px = 0; px < procShape[0]; ++px)
          {
            const Interval Gx = expandByOne(cellSplits[0][px], cx);
            const Interval Vx{
              Gx.begin,
              Gx.empty() ? Gx.end : std::min(nx, Gx.end + 1)
            };
            if (!Vx.contains(i))
              continue;

            for (int py = 0; py < procShape[1]; ++py)
            {
              const Interval Gy = expandByOne(cellSplits[1][py], cy);
              const Interval Vy{
                Gy.begin,
                Gy.empty() ? Gy.end : std::min(ny, Gy.end + 1)
              };
              if (!Vy.contains(j))
                continue;

              for (int pz = 0; pz < procShape[2]; ++pz)
              {
                const Interval Gz = expandByOne(cellSplits[2][pz], cz);
                const Interval Vz{
                  Gz.begin,
                  Gz.empty() ? Gz.end : std::min(nz, Gz.end + 1)
                };
                if (Vz.contains(k))
                  res.push_back(procCoordToRank({ px, py, pz }, procShape));
              }
            }
          }
          uniqueSort(res);
          return res;
        };

        const auto stateOfVertex = [&](size_t i, size_t j, size_t k) -> Shard::State
        {
          const int owner = ownerOfVertex(i, j, k);
          if (localVertX.contains(i) && localVertY.contains(j) && localVertZ.contains(k))
            return owner == rank ? Shard::State::Owned : Shard::State::Shared;
          return Shard::State::Ghost;
        };

        const auto holdersOfCell = [&](size_t i, size_t j, size_t k)
        {
          std::vector<int> res;
          const int owner = ownerOfCell(i, j, k);
          const std::vector<int> oc = rankToProcCoord(owner, procShape);

          for (int dx = -1; dx <= 1; ++dx)
          {
            for (int dy = -1; dy <= 1; ++dy)
            {
              for (int dz = -1; dz <= 1; ++dz)
              {
                std::vector<int> pc = oc;
                pc[0] += dx;
                pc[1] += dy;
                pc[2] += dz;
                if (pc[0] < 0 || pc[0] >= procShape[0]) continue;
                if (pc[1] < 0 || pc[1] >= procShape[1]) continue;
                if (pc[2] < 0 || pc[2] >= procShape[2]) continue;

                const Interval GX = expandByOne(cellSplits[0][pc[0]], cx);
                const Interval GY = expandByOne(cellSplits[1][pc[1]], cy);
                const Interval GZ = expandByOne(cellSplits[2][pc[2]], cz);

                if (GX.contains(i) && GY.contains(j) && GZ.contains(k))
                  res.push_back(procCoordToRank(pc, procShape));
              }
            }
          }

          uniqueSort(res);
          return res;
        };

        Shard::Builder sb;
        sb.initialize(/*dimension=*/3, /*sdim=*/3);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve(vertX.size() * vertY.size() * vertZ.size());

        for (size_t k = vertZ.begin; k < vertZ.end; ++k)
        {
          for (size_t j = vertY.begin; j < vertY.end; ++j)
          {
            for (size_t i = vertX.begin; i < vertX.end; ++i)
            {
              const Index gvid = vid3(i, j, k, nx, ny);
              const int owner = ownerOfVertex(i, j, k);
              const Shard::State state = stateOfVertex(i, j, k);

              const Index lv = sb.vertex(
                gvid,
                sp3(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k)),
                state);
              gv2lv.emplace(gvid, lv);

              if (state == Shard::State::Owned)
              {
                for (const int r : holdersOfVertex(i, j, k))
                {
                  if (r != rank)
                    sb.halo(0, lv, static_cast<Index>(r));
                }
              }
              else
              {
                sb.setOwner(0, lv, static_cast<Index>(owner));
              }
            }
          }
        }

        for (size_t k = ghostZ.begin; k < ghostZ.end; ++k)
        {
          for (size_t j = ghostY.begin; j < ghostY.end; ++j)
          {
            for (size_t i = ghostX.begin; i < ghostX.end; ++i)
            {
              const Index mid = macroId3(i, j, k, cx, cy);
              const int owner = ownerOfCell(i, j, k);
              const Shard::State state =
                (owner == rank ? Shard::State::Owned : Shard::State::Ghost);

              if (g == Polytope::Type::Hexahedron)
              {
                const Index v0 = gv2lv.at(vid3(i,     j,     k,     nx, ny));
                const Index v1 = gv2lv.at(vid3(i + 1, j,     k,     nx, ny));
                const Index v2 = gv2lv.at(vid3(i + 1, j + 1, k,     nx, ny));
                const Index v3 = gv2lv.at(vid3(i,     j + 1, k,     nx, ny));
                const Index v4 = gv2lv.at(vid3(i,     j,     k + 1, nx, ny));
                const Index v5 = gv2lv.at(vid3(i + 1, j,     k + 1, nx, ny));
                const Index v6 = gv2lv.at(vid3(i + 1, j + 1, k + 1, nx, ny));
                const Index v7 = gv2lv.at(vid3(i,     j + 1, k + 1, nx, ny));

                const IndexArray vs = makeIndexArray(v0, v1, v2, v3, v4, v5, v6, v7);
                const Index lc = sb.polytope(3, mid, Polytope::Type::Hexahedron, vs, state);

                if (state == Shard::State::Owned)
                {
                  for (const int r : holdersOfCell(i, j, k))
                  {
                    if (r != rank)
                      sb.halo(3, lc, static_cast<Index>(r));
                  }
                }
                else
                {
                  sb.setOwner(3, lc, static_cast<Index>(owner));
                }
              }
              else if (g == Polytope::Type::Wedge)
              {
                const Index v0  = gv2lv.at(vid3(i,     j,     k,     nx, ny));
                const Index v1  = gv2lv.at(vid3(i + 1, j,     k,     nx, ny));
                const Index v2  = gv2lv.at(vid3(i,     j + 1, k,     nx, ny));
                const Index v3  = gv2lv.at(vid3(i + 1, j + 1, k,     nx, ny));
                const Index v0p = gv2lv.at(vid3(i,     j,     k + 1, nx, ny));
                const Index v1p = gv2lv.at(vid3(i + 1, j,     k + 1, nx, ny));
                const Index v2p = gv2lv.at(vid3(i,     j + 1, k + 1, nx, ny));
                const Index v3p = gv2lv.at(vid3(i + 1, j + 1, k + 1, nx, ny));

                const Index lc0 = sb.polytope(
                  3, static_cast<Index>(2 * mid + 0), Polytope::Type::Wedge,
                  makeIndexArray(v0, v1, v2, v0p, v1p, v2p), state);

                const Index lc1 = sb.polytope(
                  3, static_cast<Index>(2 * mid + 1), Polytope::Type::Wedge,
                  makeIndexArray(v1, v3, v2, v1p, v3p, v2p), state);

                if (state == Shard::State::Owned)
                {
                  for (const int r : holdersOfCell(i, j, k))
                  {
                    if (r != rank)
                    {
                      sb.halo(3, lc0, static_cast<Index>(r));
                      sb.halo(3, lc1, static_cast<Index>(r));
                    }
                  }
                }
                else
                {
                  const Index ow = static_cast<Index>(owner);
                  sb.setOwner(3, lc0, ow);
                  sb.setOwner(3, lc1, ow);
                }
              }
              else // Tetrahedron
              {
                const Index v000 = gv2lv.at(vid3(i,     j,     k,     nx, ny));
                const Index v100 = gv2lv.at(vid3(i + 1, j,     k,     nx, ny));
                const Index v010 = gv2lv.at(vid3(i,     j + 1, k,     nx, ny));
                const Index v110 = gv2lv.at(vid3(i + 1, j + 1, k,     nx, ny));

                const Index v001 = gv2lv.at(vid3(i,     j,     k + 1, nx, ny));
                const Index v101 = gv2lv.at(vid3(i + 1, j,     k + 1, nx, ny));
                const Index v011 = gv2lv.at(vid3(i,     j + 1, k + 1, nx, ny));
                const Index v111 = gv2lv.at(vid3(i + 1, j + 1, k + 1, nx, ny));

                const Index lc0 = sb.polytope(
                  3, static_cast<Index>(6 * mid + 0), Polytope::Type::Tetrahedron,
                  makeIndexArray(v000, v100, v110, v111), state);
                const Index lc1 = sb.polytope(
                  3, static_cast<Index>(6 * mid + 1), Polytope::Type::Tetrahedron,
                  makeIndexArray(v000, v110, v010, v111), state);
                const Index lc2 = sb.polytope(
                  3, static_cast<Index>(6 * mid + 2), Polytope::Type::Tetrahedron,
                  makeIndexArray(v000, v010, v011, v111), state);
                const Index lc3 = sb.polytope(
                  3, static_cast<Index>(6 * mid + 3), Polytope::Type::Tetrahedron,
                  makeIndexArray(v000, v011, v001, v111), state);
                const Index lc4 = sb.polytope(
                  3, static_cast<Index>(6 * mid + 4), Polytope::Type::Tetrahedron,
                  makeIndexArray(v000, v001, v101, v111), state);
                const Index lc5 = sb.polytope(
                  3, static_cast<Index>(6 * mid + 5), Polytope::Type::Tetrahedron,
                  makeIndexArray(v000, v101, v100, v111), state);

                if (state == Shard::State::Owned)
                {
                  for (const int r : holdersOfCell(i, j, k))
                  {
                    if (r != rank)
                    {
                      sb.halo(3, lc0, static_cast<Index>(r));
                      sb.halo(3, lc1, static_cast<Index>(r));
                      sb.halo(3, lc2, static_cast<Index>(r));
                      sb.halo(3, lc3, static_cast<Index>(r));
                      sb.halo(3, lc4, static_cast<Index>(r));
                      sb.halo(3, lc5, static_cast<Index>(r));
                    }
                  }
                }
                else
                {
                  const Index ow = static_cast<Index>(owner);
                  sb.setOwner(3, lc0, ow);
                  sb.setOwner(3, lc1, ow);
                  sb.setOwner(3, lc2, ow);
                  sb.setOwner(3, lc3, ow);
                  sb.setOwner(3, lc4, ow);
                  sb.setOwner(3, lc5, ow);
                }
              }
            }
          }
        }

        return finish(sb.finalize());
      }

      default:
      {
        Alert::NamespacedException("Rodin::Geometry::Mesh<Context::MPI>::UniformGrid")
          << "Unsupported geometry type " << g << "."
          << Alert::Raise;
        return MPIMesh(context);
      }
    }
  }

}
