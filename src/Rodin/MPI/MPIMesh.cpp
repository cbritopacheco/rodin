#include "MPIMesh.h"

#ifdef RODIN_USE_MPI

namespace Rodin::Geometry
{
  MPIMesh::Builder& MPIMesh::Builder::initialize(const Context::MPI& context, Shard&& shard)
  {
    m_context = std::move(context);
    m_shard = std::move(shard);
    return *this;
  }

  // MPIMesh MPIMesh::Builder::finalize()
  // {
  //   MPIMesh mesh(m_context);
  //   mesh.m_shard = std::move(m_shard);
  //   return mesh;
  // }

  Index MPIMesh::getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId)
  {
    return 0;
  }

  MPIMesh& MPIMesh::scale(Real c)
  {
    m_shard.scale(c);
    return *this;
  }

  void MPIMesh::flush()
  {
    m_shard.flush();
  }

  bool MPIMesh::isSubMesh() const
  {
    return false;
  }

  size_t MPIMesh::getDimension() const
  {
    return m_shard.getDimension();
  }

  size_t MPIMesh::getSpaceDimension() const
  {
    return m_shard.getSpaceDimension();
  }

  const Context::MPI& MPIMesh::getContext() const
  {
    return m_context;
  }

  size_t MPIMesh::getPolytopeCount(size_t d) const
  {
    const auto& comm = m_context.getCommunicator();
    size_t count = m_shard.getPolytopeCount(d);
    boost::mpi::all_reduce(comm, count, std::plus<size_t>());
    return count;
  }

  size_t MPIMesh::getPolytopeCount(Polytope::Type g) const
  {
    const auto& comm = m_context.getCommunicator();
    size_t count = m_shard.getPolytopeCount(g);
    boost::mpi::all_reduce(comm, count, std::plus<size_t>());
    return count;
  }

  Shard& MPIMesh::getShard()
  {
    return m_shard;
  }

  const Shard& MPIMesh::getShard() const
  {
    return m_shard;
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index globalIdx) const
  {
    const auto& shard = getShard();
    std::optional<Index> localIndex = getLocalIndex(dimension, globalIdx);
    if (localIndex)
      return shard.getPolytope(dimension, *localIndex);
    else
      return PolytopeIterator(dimension, shard, BoundedIndexGenerator(0, 0));
  }

  CellIterator MPIMesh::getCell(Index globalIdx) const
  {
    const auto& shard = getShard();
    const size_t d = getDimension();
    std::optional<Index> localIndex = getLocalIndex(d, globalIdx);
    if (localIndex)
      return shard.getCell(*localIndex);
    else
      return CellIterator(shard, BoundedIndexGenerator(0, 0));
  }

  FaceIterator MPIMesh::getFace(Index globalIdx) const
  {
    const auto& shard = getShard();
    const size_t d = getDimension() - 1;
    std::optional<Index> localIndex = getLocalIndex(d, globalIdx);
    if (localIndex)
      return shard.getFace(*localIndex);
    else
      return FaceIterator(shard, BoundedIndexGenerator(0, 0));
  }

  VertexIterator MPIMesh::getVertex(Index globalIdx) const
  {
    const auto& shard = getShard();
    const size_t d = getDimension() - 1;
    std::optional<Index> localIndex = getLocalIndex(d, globalIdx);
    if (localIndex)
      return shard.getVertex(*localIndex);
    else
      return VertexIterator(shard, BoundedIndexGenerator(0, 0));
  }

  std::optional<Index> MPIMesh::getLocalIndex(size_t dimension, Index globalIdx) const
  {
    const auto& shard = getShard();
    const auto& polytopeMap = shard.getPolytopeMap(dimension);
    const auto find = polytopeMap.right.find(globalIdx);
    if (find == polytopeMap.right.end())
      return std::nullopt;
    else
      return find->get_left();
  }

  bool MPIMesh::isInterface(Index faceIdx) const
  {
    const auto& comm = m_context.getCommunicator();
    const size_t D = getDimension();
    const auto& shard = getShard();
    const auto& conn = shard.getConnectivity();
    const auto localIndex = getLocalIndex(D - 1, faceIdx);
    bool result = false;
    if (localIndex)
    {
      RODIN_GEOMETRY_REQUIRE_INCIDENCE(shard, D - 1, D);
      assert(conn.getIncidence(D - 1, D).size());
      const auto& incidence = conn.getIncidence({D - 1, D}, *localIndex);
      assert(incidence.size() > 0);
      result = (incidence.size() > 1);
    }
    boost::mpi::all_reduce(comm, result, std::logical_or<bool>());
    return result;
  }

  bool MPIMesh::isBoundary(Index faceIdx) const
  {
    const auto& comm = m_context.getCommunicator();
    const size_t D = getDimension();
    const auto& shard = getShard();
    const auto& conn = shard.getConnectivity();
    const auto localIndex = getLocalIndex(D - 1, faceIdx);
    bool result = false;
    if (localIndex)
    {
      RODIN_GEOMETRY_REQUIRE_INCIDENCE(shard, D - 1, D);
      assert(conn.getIncidence(D - 1, D).size());
      const auto& incidence = conn.getIncidence({D - 1, D}, *localIndex);
      assert(incidence.size() > 0);
      result = (incidence.size() == 1);
    }
    boost::mpi::all_reduce(comm, result, std::logical_or<bool>());
    return result;
  }
}

#endif
