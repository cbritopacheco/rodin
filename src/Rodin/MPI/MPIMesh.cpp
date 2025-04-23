#include <boost/serialization/optional.hpp>

#include "MPIMesh.h"

namespace Rodin::Geometry
{
  MPIMesh::Builder& MPIMesh::Builder::initialize(const Context::MPI& context, Shard&& shard)
  {
    m_context = context;
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
    getShard().scale(c);
    return *this;
  }

  void MPIMesh::flush()
  {
    getShard().flush();
  }

  bool MPIMesh::isSubMesh() const
  {
    return false;
  }

  size_t MPIMesh::getDimension() const
  {
    return getShard().getDimension();
  }

  size_t MPIMesh::getSpaceDimension() const
  {
    return getShard().getSpaceDimension();
  }

  const Context::MPI& MPIMesh::getContext() const
  {
    return m_context;
  }

  size_t MPIMesh::getPolytopeCount(size_t d) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    size_t local = shard.getPolytopeCount(d) - shard.getGhosts()[d].size();
    boost::mpi::all_reduce(comm, local, std::plus<size_t>());
    return local;
  }

  size_t MPIMesh::getPolytopeCount(Polytope::Type g) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    const size_t d = Polytope::getGeometryDimension(g);
    size_t local = 0;
    for (auto it = shard.getPolytope(d); it; ++it)
    {
      if (!shard.isGhost(d, it->getIndex()))
      {
        if (it->getGeometry() == g)
          local++;
      }
    }
    boost::mpi::all_reduce(comm, local, std::plus<size_t>());
    return local;
  }

  Shard& MPIMesh::getShard()
  {
    return m_shard;
  }

  const Shard& MPIMesh::getShard() const
  {
    return m_shard;
  }

  std::optional<Index> MPIMesh::getLocalIndex(size_t dimension, Index globalIdx) const
  {
    const auto& map = getShard().getPolytopeMap(dimension).right;
    auto it = map.find(globalIdx);
    if (it == map.end())
      return std::nullopt;
    return it->get_left();
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index globalIdx) const
  {
    auto local = getLocalIndex(dimension, globalIdx);
    if (local)
      return getShard().getPolytope(dimension, *local);
    return PolytopeIterator(dimension, getShard(), BoundedIndexGenerator(0, 0));
  }

  CellIterator MPIMesh::getCell(Index globalIdx) const
  {
    const size_t d = getDimension();
    auto local = getLocalIndex(d, globalIdx);
    if (local)
      return getShard().getCell(*local);
    return CellIterator(getShard(), BoundedIndexGenerator(0, 0));
  }

  FaceIterator MPIMesh::getFace(Index globalIdx) const
  {
    const size_t d = getDimension() - 1;
    auto local = getLocalIndex(d, globalIdx);
    if (local)
      return getShard().getFace(*local);
    return FaceIterator(getShard(), BoundedIndexGenerator(0, 0));
  }

  VertexIterator MPIMesh::getVertex(Index globalIdx) const
  {
    auto local = getLocalIndex(0, globalIdx);
    if (local)
      return getShard().getVertex(*local);
    return VertexIterator(getShard(), BoundedIndexGenerator(0, 0));
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension) const
  {
    return getShard().getPolytope(dimension);
  }

  CellIterator MPIMesh::getCell() const
  {
    return getShard().getCell();
  }

  FaceIterator MPIMesh::getFace() const
  {
    return getShard().getFace();
  }

  VertexIterator MPIMesh::getVertex() const
  {
    return getShard().getVertex();
  }

  FaceIterator MPIMesh::getBoundary() const
  {
    return getShard().getBoundary();
  }

  FaceIterator MPIMesh::getInterface() const
  {
    return getShard().getInterface();
  }

  bool MPIMesh::isInterface(Index faceIdx) const
  {
    const auto& comm = m_context.getCommunicator();
    const size_t D = getDimension();
    auto local = getLocalIndex(D - 1, faceIdx);
    bool res = false;
    if (local)
    {
      RODIN_GEOMETRY_REQUIRE_INCIDENCE(getShard(), D - 1, D);
      const auto& inc = getShard().getConnectivity()
                                  .getIncidence({D - 1, D}, *local);
      res = inc.size() > 1;
    }
    boost::mpi::all_reduce(comm, res, std::logical_or<bool>());
    return res;
  }

  bool MPIMesh::isBoundary(Index faceIdx) const
  {
    const auto& comm = m_context.getCommunicator();
    const size_t D = getDimension();
    auto local = getLocalIndex(D - 1, faceIdx);
    bool res = false;
    if (local)
    {
      RODIN_GEOMETRY_REQUIRE_INCIDENCE(getShard(), D - 1, D);
      const auto& inc = getShard().getConnectivity()
                                  .getIncidence({D - 1, D}, *local);
      res = inc.size() == 1;
    }
    boost::mpi::all_reduce(comm, res, std::logical_or<bool>());
    return res;
  }

  Real MPIMesh::getMeasure(size_t d) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getMeasure(d);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getMeasure(size_t d, Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getMeasure(d, attr);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getMeasure(size_t d, const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getMeasure(d, attrs);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getVolume() const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getVolume();
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getVolume(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getVolume(attr);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getVolume(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getVolume(attrs);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getArea() const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getArea();
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getArea(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getArea(attr);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getArea(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getArea(attrs);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getPerimeter() const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getPerimeter();
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getPerimeter(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getPerimeter(attr);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  Real MPIMesh::getPerimeter(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    Real local = getShard().getPerimeter(attrs);
    boost::mpi::all_reduce(comm, local, std::plus<Real>());
    return local;
  }

  const PolytopeTransformation& MPIMesh::getPolytopeTransformation(size_t dimension, Index globalIdx) const
  {
    auto local = getLocalIndex(dimension, globalIdx);
    if (local)
    {
      const auto& comm = m_context.getCommunicator();
    }
    return getShard().getPolytopeTransformation(dimension, *local);
  }

  Polytope::Type MPIMesh::getGeometry(size_t dimension, Index globalIdx) const
  {
    auto local = getLocalIndex(dimension, globalIdx);
    boost::optional<Polytope::Type> res;
    if (local)
    {
      const auto& comm = m_context.getCommunicator();
      res = getShard().getGeometry(dimension, *local);
      boost::mpi::broadcast(comm, res, comm.rank());
    }
    assert(res);
    return *res;
  }

  Attribute MPIMesh::getAttribute(size_t dimension, Index globalIdx) const
  {
    auto local = getLocalIndex(dimension, globalIdx);
    boost::optional<Attribute> res;
    if (local)
    {
      const auto& comm = m_context.getCommunicator();
      res = getShard().getAttribute(dimension, *local);
      boost::mpi::broadcast(comm, res, comm.rank());
    }
    assert(res);
    return *res;
  }

  MPIMesh& MPIMesh::setAttribute(const std::pair<size_t, Index>& p, Attribute attr)
  {
    auto local = getLocalIndex(p.first, p.second);
    if (local)
      getShard().setAttribute({ p.first, *local }, attr);
    return *this;
  }

  void MPIMesh::getVertexCoordinates(Math::SpatialVector<Real>& out, Index globalIdx) const
  {
    auto local = getLocalIndex(0, globalIdx);
    if (local)
    {
      const auto& comm = m_context.getCommunicator();
      getShard().getVertexCoordinates(out, *local);
      boost::mpi::broadcast(comm, out, comm.rank());
    }
    assert(out.size() == getSpaceDimension());
  }

  Real MPIMesh::getVertexCoordinates(Index globalIdx, size_t i) const
  {
    auto local = getLocalIndex(0, globalIdx);
    boost::optional<Real> res;
    if (local)
    {
      const auto& comm = m_context.getCommunicator();
      res = getShard().getVertexCoordinates(*local, i);
      boost::mpi::broadcast(comm, res, comm.rank());
    }
    assert(res);
    return *res;
  }

  MPIMesh& MPIMesh::setVertexCoordinates(Index globalIdx, Real s, size_t i)
  {
    auto local = getLocalIndex(0, globalIdx);
    if (local)
      getShard().setVertexCoordinates(*local, s, i);
    return *this;
  }

  MPIMesh& MPIMesh::setVertexCoordinates(Index globalIdx, const Math::SpatialVector<Real>& coords)
  {
    auto local = getLocalIndex(0, globalIdx);
    if (local)
      getShard().setVertexCoordinates(*local, coords);
    return *this;
  }

  MPIMesh& MPIMesh::setPolytopeTransformation(
      const std::pair<size_t, Index> p, PolytopeTransformation* trans)
  {
    auto local = getLocalIndex(p.first, p.second);
    if (local)
      getShard().setPolytopeTransformation({ p.first, *local }, trans);
    return *this;
  }
}

