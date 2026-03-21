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

  MPIMesh MPIMesh::UniformGrid(
    const Context::MPI& context, Polytope::Type g, const Array<size_t>& shape)
  {
    const auto& comm = context.getCommunicator();
    const int rank = comm.rank();
    const int size = comm.size();

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

    const auto splitRange = [](size_t n, int r, int p)
    {
      struct Interval
      {
        size_t begin;
        size_t end;
      };

      const size_t q = n / static_cast<size_t>(p);
      const size_t rem = n % static_cast<size_t>(p);
      const size_t begin =
        static_cast<size_t>(r) * q + std::min(static_cast<size_t>(r), rem);
      const size_t count = q + (static_cast<size_t>(r) < rem ? 1 : 0);
      return Interval{ begin, begin + count };
    };

    const auto buildColumnOwner = [&](size_t cx)
    {
      std::vector<int> owner(cx, -1);
      for (int r = 0; r < size; ++r)
      {
        const auto I = splitRange(cx, r, size);
        for (size_t i = I.begin; i < I.end; ++i)
          owner[i] = r;
      }
      return owner;
    };

    const auto ownerOfVertexX = [](size_t i, const std::vector<int>& cellOwner)
    {
      int owner = std::numeric_limits<int>::max();

      if (i > 0 && i - 1 < cellOwner.size())
        owner = std::min(owner, cellOwner[i - 1]);

      if (i < cellOwner.size())
        owner = std::min(owner, cellOwner[i]);

      if (owner == std::numeric_limits<int>::max())
        owner = 0;

      return owner;
    };

    const auto sharersOfVertexX = [](size_t i, const std::vector<int>& cellOwner)
    {
      std::vector<int> res;

      const auto add = [&](int r)
      {
        if (r < 0)
          return;
        if (std::find(res.begin(), res.end(), r) == res.end())
          res.push_back(r);
      };

      if (i > 0 && i - 1 < cellOwner.size())
        add(cellOwner[i - 1]);

      if (i < cellOwner.size())
        add(cellOwner[i]);

      std::sort(res.begin(), res.end());
      return res;
    };

    const auto sharersOfCellColumn = [](size_t i, const std::vector<int>& cellOwner)
    {
      std::vector<int> res;
      const int owner = cellOwner[i];

      if (i > 0 && cellOwner[i - 1] != owner)
        res.push_back(cellOwner[i - 1]);

      if (i + 1 < cellOwner.size() && cellOwner[i + 1] != owner)
      {
        if (std::find(res.begin(), res.end(), cellOwner[i + 1]) == res.end())
          res.push_back(cellOwner[i + 1]);
      }

      std::sort(res.begin(), res.end());
      return res;
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

        const Shard::Flags flags =
          (rank == 0 ? Shard::Flags::Owned : Shard::Flags::Ghost);

        const Index lv = sb.vertex(/*globalIdx=*/0, sp0(), flags);

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

        const size_t cx = nx - 1;
        const auto xr = splitRange(cx, rank, size);
        const std::vector<int> cellOwner = buildColumnOwner(cx);

        size_t gx0 = xr.begin;
        size_t gx1 = xr.end;
        if (cx > 0)
        {
          if (gx0 > 0) gx0--;
          if (gx1 < cx) gx1++;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/1, /*sdim=*/1);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve(gx1 - gx0 + 1);

        for (size_t i = gx0; i <= gx1; ++i)
        {
          const Index gvid = vid1(i);
          const int owner = ownerOfVertexX(i, cellOwner);
          const Shard::Flags flags =
            (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

          const Index lv = sb.vertex(gvid, sp1(static_cast<Real>(i)), flags);
          gv2lv.emplace(gvid, lv);

          if (owner == rank)
          {
            for (const int r : sharersOfVertexX(i, cellOwner))
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

        for (size_t i = gx0; i < gx1; ++i)
        {
          const int owner = cellOwner[i];
          const Shard::Flags flags =
            (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

          const Index gcid = static_cast<Index>(i);
          const IndexArray vs = makeIndexArray(
            gv2lv.at(vid1(i)),
            gv2lv.at(vid1(i + 1))
          );

          const Index lc = sb.polytope(1, gcid, Polytope::Type::Segment, vs, flags);

          if (owner == rank)
          {
            for (const int r : sharersOfCellColumn(i, cellOwner))
              sb.halo(1, lc, static_cast<Index>(r));
          }
          else
          {
            sb.setOwner(1, lc, static_cast<Index>(owner));
          }
        }

        return finish(sb.finalize());
      }

      case Polytope::Type::Triangle:
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
            << "Triangle uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;
        const auto xr = splitRange(cx, rank, size);
        const std::vector<int> cellOwner = buildColumnOwner(cx);

        size_t gx0 = xr.begin;
        size_t gx1 = xr.end;
        if (cx > 0)
        {
          if (gx0 > 0) gx0--;
          if (gx1 < cx) gx1++;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/2, /*sdim=*/2);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve((gx1 - gx0 + 1) * ny);

        for (size_t j = 0; j < ny; ++j)
        {
          for (size_t i = gx0; i <= gx1; ++i)
          {
            const Index gvid = vid2(i, j, nx);
            const int owner = ownerOfVertexX(i, cellOwner);
            const Shard::Flags flags =
              (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

            const Index lv = sb.vertex(
              gvid,
              sp2(static_cast<Real>(i), static_cast<Real>(j)),
              flags);
            gv2lv.emplace(gvid, lv);

            if (owner == rank)
            {
              for (const int r : sharersOfVertexX(i, cellOwner))
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

        for (size_t j = 0; j < cy; ++j)
        {
          for (size_t i = gx0; i < gx1; ++i)
          {
            const Index mid = macroId2(i, j, cx);
            const int owner = cellOwner[i];
            const Shard::Flags flags =
              (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

            const Index v00 = gv2lv.at(vid2(i,     j,     nx));
            const Index v10 = gv2lv.at(vid2(i + 1, j,     nx));
            const Index v01 = gv2lv.at(vid2(i,     j + 1, nx));
            const Index v11 = gv2lv.at(vid2(i + 1, j + 1, nx));

            {
              const Index gcid = static_cast<Index>(2 * mid);
              const IndexArray vs = makeIndexArray(v00, v10, v01);
              const Index lc = sb.polytope(2, gcid, Polytope::Type::Triangle, vs, flags);

              if (owner == rank)
              {
                for (const int r : sharersOfCellColumn(i, cellOwner))
                  sb.halo(2, lc, static_cast<Index>(r));
              }
              else
              {
                sb.setOwner(2, lc, static_cast<Index>(owner));
              }
            }

            {
              const Index gcid = static_cast<Index>(2 * mid + 1);
              const IndexArray vs = makeIndexArray(v10, v11, v01);
              const Index lc = sb.polytope(2, gcid, Polytope::Type::Triangle, vs, flags);

              if (owner == rank)
              {
                for (const int r : sharersOfCellColumn(i, cellOwner))
                  sb.halo(2, lc, static_cast<Index>(r));
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
            << "Quadrilateral uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;
        const auto xr = splitRange(cx, rank, size);
        const std::vector<int> cellOwner = buildColumnOwner(cx);

        size_t gx0 = xr.begin;
        size_t gx1 = xr.end;
        if (cx > 0)
        {
          if (gx0 > 0) gx0--;
          if (gx1 < cx) gx1++;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/2, /*sdim=*/2);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve((gx1 - gx0 + 1) * ny);

        for (size_t j = 0; j < ny; ++j)
        {
          for (size_t i = gx0; i <= gx1; ++i)
          {
            const Index gvid = vid2(i, j, nx);
            const int owner = ownerOfVertexX(i, cellOwner);
            const Shard::Flags flags =
              (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

            const Index lv = sb.vertex(
              gvid,
              sp2(static_cast<Real>(i), static_cast<Real>(j)),
              flags);
            gv2lv.emplace(gvid, lv);

            if (owner == rank)
            {
              for (const int r : sharersOfVertexX(i, cellOwner))
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

        for (size_t j = 0; j < cy; ++j)
        {
          for (size_t i = gx0; i < gx1; ++i)
          {
            const Index gcid = macroId2(i, j, cx);
            const int owner = cellOwner[i];
            const Shard::Flags flags =
              (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

            const IndexArray vs = makeIndexArray(
              gv2lv.at(vid2(i,     j,     nx)),
              gv2lv.at(vid2(i + 1, j,     nx)),
              gv2lv.at(vid2(i + 1, j + 1, nx)),
              gv2lv.at(vid2(i,     j + 1, nx))
            );

            const Index lc = sb.polytope(2, gcid, Polytope::Type::Quadrilateral, vs, flags);

            if (owner == rank)
            {
              for (const int r : sharersOfCellColumn(i, cellOwner))
                sb.halo(2, lc, static_cast<Index>(r));
            }
            else
            {
              sb.setOwner(2, lc, static_cast<Index>(owner));
            }
          }
        }

        return finish(sb.finalize());
      }

      case Polytope::Type::Tetrahedron:
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
            << "Tetrahedron uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;
        const size_t cz = nz - 1;
        const auto xr = splitRange(cx, rank, size);
        const std::vector<int> cellOwner = buildColumnOwner(cx);

        size_t gx0 = xr.begin;
        size_t gx1 = xr.end;
        if (cx > 0)
        {
          if (gx0 > 0) gx0--;
          if (gx1 < cx) gx1++;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/3, /*sdim=*/3);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve((gx1 - gx0 + 1) * ny * nz);

        for (size_t k = 0; k < nz; ++k)
        {
          for (size_t j = 0; j < ny; ++j)
          {
            for (size_t i = gx0; i <= gx1; ++i)
            {
              const Index gvid = vid3(i, j, k, nx, ny);
              const int owner = ownerOfVertexX(i, cellOwner);
              const Shard::Flags flags =
                (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

              const Index lv = sb.vertex(
                gvid,
                sp3(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k)),
                flags);
              gv2lv.emplace(gvid, lv);

              if (owner == rank)
              {
                for (const int r : sharersOfVertexX(i, cellOwner))
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

        for (size_t k = 0; k < cz; ++k)
        {
          for (size_t j = 0; j < cy; ++j)
          {
            for (size_t i = gx0; i < gx1; ++i)
            {
              const Index mid = macroId3(i, j, k, cx, cy);
              const int owner = cellOwner[i];
              const Shard::Flags flags =
                (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

              const Index v000 = gv2lv.at(vid3(i,     j,     k,     nx, ny));
              const Index v100 = gv2lv.at(vid3(i + 1, j,     k,     nx, ny));
              const Index v010 = gv2lv.at(vid3(i,     j + 1, k,     nx, ny));
              const Index v110 = gv2lv.at(vid3(i + 1, j + 1, k,     nx, ny));

              const Index v001 = gv2lv.at(vid3(i,     j,     k + 1, nx, ny));
              const Index v101 = gv2lv.at(vid3(i + 1, j,     k + 1, nx, ny));
              const Index v011 = gv2lv.at(vid3(i,     j + 1, k + 1, nx, ny));
              const Index v111 = gv2lv.at(vid3(i + 1, j + 1, k + 1, nx, ny));

              const IndexArray tet0 = makeIndexArray(v000, v100, v110, v111);
              const IndexArray tet1 = makeIndexArray(v000, v110, v010, v111);
              const IndexArray tet2 = makeIndexArray(v000, v010, v011, v111);
              const IndexArray tet3 = makeIndexArray(v000, v011, v001, v111);
              const IndexArray tet4 = makeIndexArray(v000, v001, v101, v111);
              const IndexArray tet5 = makeIndexArray(v000, v101, v100, v111);

              const Index lc0 = sb.polytope(3, static_cast<Index>(6 * mid + 0), Polytope::Type::Tetrahedron, tet0, flags);
              const Index lc1 = sb.polytope(3, static_cast<Index>(6 * mid + 1), Polytope::Type::Tetrahedron, tet1, flags);
              const Index lc2 = sb.polytope(3, static_cast<Index>(6 * mid + 2), Polytope::Type::Tetrahedron, tet2, flags);
              const Index lc3 = sb.polytope(3, static_cast<Index>(6 * mid + 3), Polytope::Type::Tetrahedron, tet3, flags);
              const Index lc4 = sb.polytope(3, static_cast<Index>(6 * mid + 4), Polytope::Type::Tetrahedron, tet4, flags);
              const Index lc5 = sb.polytope(3, static_cast<Index>(6 * mid + 5), Polytope::Type::Tetrahedron, tet5, flags);

              if (owner == rank)
              {
                const auto rs = sharersOfCellColumn(i, cellOwner);
                for (const int r : rs)
                {
                  sb.halo(3, lc0, static_cast<Index>(r));
                  sb.halo(3, lc1, static_cast<Index>(r));
                  sb.halo(3, lc2, static_cast<Index>(r));
                  sb.halo(3, lc3, static_cast<Index>(r));
                  sb.halo(3, lc4, static_cast<Index>(r));
                  sb.halo(3, lc5, static_cast<Index>(r));
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

        return finish(sb.finalize());
      }

      case Polytope::Type::Hexahedron:
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
            << "Hexahedron uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;
        const size_t cz = nz - 1;
        const auto xr = splitRange(cx, rank, size);
        const std::vector<int> cellOwner = buildColumnOwner(cx);

        size_t gx0 = xr.begin;
        size_t gx1 = xr.end;
        if (cx > 0)
        {
          if (gx0 > 0) gx0--;
          if (gx1 < cx) gx1++;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/3, /*sdim=*/3);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve((gx1 - gx0 + 1) * ny * nz);

        for (size_t k = 0; k < nz; ++k)
        {
          for (size_t j = 0; j < ny; ++j)
          {
            for (size_t i = gx0; i <= gx1; ++i)
            {
              const Index gvid = vid3(i, j, k, nx, ny);
              const int owner = ownerOfVertexX(i, cellOwner);
              const Shard::Flags flags =
                (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

              const Index lv = sb.vertex(
                gvid,
                sp3(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k)),
                flags);
              gv2lv.emplace(gvid, lv);

              if (owner == rank)
              {
                for (const int r : sharersOfVertexX(i, cellOwner))
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

        for (size_t k = 0; k < cz; ++k)
        {
          for (size_t j = 0; j < cy; ++j)
          {
            for (size_t i = gx0; i < gx1; ++i)
            {
              const Index gcid = macroId3(i, j, k, cx, cy);
              const int owner = cellOwner[i];
              const Shard::Flags flags =
                (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

              const Index v0 = gv2lv.at(vid3(i,     j,     k,     nx, ny));
              const Index v1 = gv2lv.at(vid3(i + 1, j,     k,     nx, ny));
              const Index v2 = gv2lv.at(vid3(i + 1, j + 1, k,     nx, ny));
              const Index v3 = gv2lv.at(vid3(i,     j + 1, k,     nx, ny));
              const Index v4 = gv2lv.at(vid3(i,     j,     k + 1, nx, ny));
              const Index v5 = gv2lv.at(vid3(i + 1, j,     k + 1, nx, ny));
              const Index v6 = gv2lv.at(vid3(i + 1, j + 1, k + 1, nx, ny));
              const Index v7 = gv2lv.at(vid3(i,     j + 1, k + 1, nx, ny));

              const IndexArray vs = makeIndexArray(v0, v1, v2, v3, v4, v5, v6, v7);
              const Index lc = sb.polytope(3, gcid, Polytope::Type::Hexahedron, vs, flags);

              if (owner == rank)
              {
                for (const int r : sharersOfCellColumn(i, cellOwner))
                  sb.halo(3, lc, static_cast<Index>(r));
              }
              else
              {
                sb.setOwner(3, lc, static_cast<Index>(owner));
              }
            }
          }
        }

        return finish(sb.finalize());
      }

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
            << "Wedge uniform grid requires at least 2 vertices per direction."
            << Alert::Raise;
        }

        const size_t cx = nx - 1;
        const size_t cy = ny - 1;
        const size_t cz = nz - 1;
        const auto xr = splitRange(cx, rank, size);
        const std::vector<int> cellOwner = buildColumnOwner(cx);

        size_t gx0 = xr.begin;
        size_t gx1 = xr.end;
        if (cx > 0)
        {
          if (gx0 > 0) gx0--;
          if (gx1 < cx) gx1++;
        }

        Shard::Builder sb;
        sb.initialize(/*dimension=*/3, /*sdim=*/3);

        UnorderedMap<Index, Index> gv2lv;
        gv2lv.reserve((gx1 - gx0 + 1) * ny * nz);

        for (size_t k = 0; k < nz; ++k)
        {
          for (size_t j = 0; j < ny; ++j)
          {
            for (size_t i = gx0; i <= gx1; ++i)
            {
              const Index gvid = vid3(i, j, k, nx, ny);
              const int owner = ownerOfVertexX(i, cellOwner);
              const Shard::Flags flags =
                (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

              const Index lv = sb.vertex(
                gvid,
                sp3(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k)),
                flags);
              gv2lv.emplace(gvid, lv);

              if (owner == rank)
              {
                for (const int r : sharersOfVertexX(i, cellOwner))
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

        for (size_t k = 0; k < cz; ++k)
        {
          for (size_t j = 0; j < cy; ++j)
          {
            for (size_t i = gx0; i < gx1; ++i)
            {
              const Index mid = macroId3(i, j, k, cx, cy);
              const int owner = cellOwner[i];
              const Shard::Flags flags =
                (owner == rank ? Shard::Flags::Owned : Shard::Flags::Ghost);

              const Index v0  = gv2lv.at(vid3(i,     j,     k,     nx, ny));
              const Index v1  = gv2lv.at(vid3(i + 1, j,     k,     nx, ny));
              const Index v2  = gv2lv.at(vid3(i,     j + 1, k,     nx, ny));
              const Index v3  = gv2lv.at(vid3(i + 1, j + 1, k,     nx, ny));
              const Index v0p = gv2lv.at(vid3(i,     j,     k + 1, nx, ny));
              const Index v1p = gv2lv.at(vid3(i + 1, j,     k + 1, nx, ny));
              const Index v2p = gv2lv.at(vid3(i,     j + 1, k + 1, nx, ny));
              const Index v3p = gv2lv.at(vid3(i + 1, j + 1, k + 1, nx, ny));

              const IndexArray w0 = makeIndexArray(v0, v1, v2, v0p, v1p, v2p);
              const IndexArray w1 = makeIndexArray(v1, v3, v2, v1p, v3p, v2p);

              const Index lc0 = sb.polytope(3, static_cast<Index>(2 * mid + 0), Polytope::Type::Wedge, w0, flags);
              const Index lc1 = sb.polytope(3, static_cast<Index>(2 * mid + 1), Polytope::Type::Wedge, w1, flags);

              if (owner == rank)
              {
                const auto rs = sharersOfCellColumn(i, cellOwner);
                for (const int r : rs)
                {
                  sb.halo(3, lc0, static_cast<Index>(r));
                  sb.halo(3, lc1, static_cast<Index>(r));
                }
              }
              else
              {
                const Index ow = static_cast<Index>(owner);
                sb.setOwner(3, lc0, ow);
                sb.setOwner(3, lc1, ow);
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

