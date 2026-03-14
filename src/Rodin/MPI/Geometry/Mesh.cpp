#include <boost/dynamic_bitset.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

#include "Rodin/Geometry/Mesh.h"
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
    return CellIterator(shard, BoundedIndexGenerator(localIdx, shard.getCellCount()));
  }

  FaceIterator MPIMesh::getFace(Index localIdx) const
  {
    const auto& shard = this->getShard();
    return FaceIterator(shard, BoundedIndexGenerator(localIdx, shard.getFaceCount()));
  }

  VertexIterator MPIMesh::getVertex(Index localIdx) const
  {
    const auto& shard = this->getShard();
    return VertexIterator(shard, BoundedIndexGenerator(localIdx, shard.getVertexCount()));
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index localIdx) const
  {
    const auto& shard = this->getShard();
    return PolytopeIterator(
        dimension, shard, BoundedIndexGenerator(localIdx, shard.getPolytopeCount(dimension)));
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
    return FaceIterator(shard, VectorIndexGenerator(std::move(indices)));
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
    return FaceIterator(shard, VectorIndexGenerator(std::move(indices)));
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

    const auto& comm = m_context.getCommunicator();
    const int rank   = comm.rank();

    const size_t D = shard.getDimension();
    assert(d <= D);

    // Vertices and cells are assumed to already carry distributed ids/ownership.
    if (d == 0 || d == D)
      return *this;

    // Make sure the two relations needed by reconciliation exist.
    // Restrict mode avoids discovering anything new here.
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(shard, D, d);
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(shard, d, 0);

    const auto& D2d = conn.getIncidence(D, d);
    const auto& d20 = conn.getIncidence(d, 0);

    const size_t nd = shard.getPolytopeCount(d);
    const size_t nc = shard.getPolytopeCount(D);

    if (nd == 0)
      return *this;

    // --------------------------------------------------------------------------
    // Neighbor stencil from cell ownership metadata.
    // --------------------------------------------------------------------------
    std::vector<int> neighbors;
    {
      UnorderedSet<int> nbrs;

      const auto& cellHalo  = shard.getHalo(D);   // owned cell -> remote ranks
      const auto& cellOwner = shard.getOwner(D);  // ghost cell -> owner rank

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
    }

    // --------------------------------------------------------------------------
    // Build a canonical key for each local d-entity from distributed vertex ids.
    //
    // Since reconcile(d) is done for a fixed d, the distributed vertex key is
    // enough for the currently supported polytope families.
    // --------------------------------------------------------------------------
    std::vector<Polytope::Key> keys(nd);
    UnorderedMap<Polytope::Key, Index, Polytope::Key::SymmetricHash, Polytope::Key::SymmetricEquality> key2local;
    key2local.reserve(nd);

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      const auto& lv = d20[i];
      Polytope::Key key(static_cast<std::uint8_t>(lv.size()));

      for (std::uint8_t j = 0; j < static_cast<std::uint8_t>(lv.size()); ++j)
      {
        const Index localVertex = lv[j];
        key[j] = shard.getPolytopeMap(0).left.at(localVertex);
      }

      keys[i] = key;
      key2local.emplace(keys[i], i);
    }

    // --------------------------------------------------------------------------
    // Local classification from incident cells.
    //
    // ownerCandidates[i] is the set of ranks that may own entity i:
    //   - self, if i touches at least one owned cell;
    //   - owner(ghost-cell), for each incident ghost cell.
    //
    // This is enough to decide ownership without first exchanging the keys.
    // --------------------------------------------------------------------------
    std::vector<uint8_t> hasOwnedIncident(nd, 0);
    std::vector<int> ownerRank(nd, std::numeric_limits<int>::max());
    std::vector<UnorderedSet<int>> participants(nd);

    for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
    {
      const bool owned = shard.isOwned(D, cell);
      const bool ghost = shard.isGhost(D, cell);

      int ghostOwner = std::numeric_limits<int>::max();
      if (ghost)
      {
        const auto it = shard.getOwner(D).find(cell);
        assert(it != shard.getOwner(D).end());
        ghostOwner = static_cast<int>(it->second);
      }

      for (const Index e : D2d[cell])
      {
        if (owned)
        {
          hasOwnedIncident[e] = 1;
          participants[e].insert(rank);
        }

        if (ghost)
        {
          participants[e].insert(ghostOwner);
        }
      }
    }

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
    {
      assert(!participants[i].empty());

      int o = std::numeric_limits<int>::max();
      for (const int r : participants[i])
        o = std::min(o, r);

      ownerRank[i] = o;
    }

    // --------------------------------------------------------------------------
    // Assign distributed ids to locally owned entities.
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
    // Owner -> ghost holders: send back distributed ids.
    //
    // Deadlock-free because every rank sends exactly one message (possibly empty)
    // to every neighbor and receives exactly one from every neighbor.
    // --------------------------------------------------------------------------
    const int tagIds = 2000 + static_cast<int>(d);

    {
      std::vector<std::vector<std::pair<Polytope::Key, Index>>> sendbuf(neighbors.size());

      for (Index i = 0; i < static_cast<Index>(nd); ++i)
      {
        if (ownerRank[i] != rank)
          continue;

        for (size_t k = 0; k < neighbors.size(); ++k)
        {
          const int nbr = neighbors[k];
          if (participants[i].find(nbr) != participants[i].end())
            sendbuf[k].push_back({ keys[i], distId[i] });
        }
      }

      std::vector<boost::mpi::request> reqs;
      reqs.reserve(neighbors.size());

      for (size_t k = 0; k < neighbors.size(); ++k)
        reqs.push_back(comm.isend(neighbors[k], tagIds, sendbuf[k]));

      for (size_t k = 0; k < neighbors.size(); ++k)
      {
        std::vector<std::pair<Polytope::Key, Index>> recvbuf;
        comm.recv(neighbors[k], tagIds, recvbuf);

        for (const auto& [key, gid] : recvbuf)
        {
          const auto it = key2local.find(key);
          if (it == key2local.end())
            continue;

          const Index i = it->second;
          if (ownerRank[i] != rank)
            distId[i] = gid;
        }
      }

      boost::mpi::wait_all(reqs.begin(), reqs.end());
    }

    for (Index i = 0; i < static_cast<Index>(nd); ++i)
      assert(distId[i] != std::numeric_limits<Index>::max());

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

      auto& flags = shard.getFlags(d);
      flags.assign(nd, Shard::Flags::None);

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
          for (const int r : participants[i])
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

