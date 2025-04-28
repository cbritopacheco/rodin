#include <boost/serialization/optional.hpp>

#include "MPIMesh.h"

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
    size_t res = shard.getPolytopeCount(d) - shard.getGhosts()[d].size();
    boost::mpi::all_reduce(comm, res, std::plus<size_t>());
    return res;
  }

  size_t MPIMesh::getPolytopeCount(Polytope::Type g) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    const size_t d = Polytope::getGeometryDimension(g);
    size_t res = 0;
    for (auto it = shard.getPolytope(d); it; ++it)
    {
      if (!shard.isGhost(d, it->getIndex()))
      {
        if (it->getGeometry() == g)
          res++;
      }
    }
    boost::mpi::all_reduce(comm, res, std::plus<size_t>());
    return res;
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

  CellIterator MPIMesh::getCell() const
  {
    return getCell(0);
  }

  FaceIterator MPIMesh::getFace() const
  {
    return getFace(0);
  }

  VertexIterator MPIMesh::getVertex() const
  {
    return getVertex(0);
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension) const
  {
    return getPolytope(dimension, 0);
  }

  CellIterator MPIMesh::getCell(Index globalIdx) const
  {
    return CellIterator(*this, BoundedIndexGenerator(globalIdx, getCellCount()));
  }

  FaceIterator MPIMesh::getFace(Index globalIdx) const
  {
    return FaceIterator(*this, BoundedIndexGenerator(globalIdx, getFaceCount()));
  }

  VertexIterator MPIMesh::getVertex(Index globalIdx) const
  {
    return VertexIterator(*this, BoundedIndexGenerator(globalIdx, getVertexCount()));
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index globalIdx) const
  {
    return PolytopeIterator(dimension, *this, BoundedIndexGenerator(globalIdx, getPolytopeCount(dimension)));
  }

  FaceIterator MPIMesh::getBoundary() const
  {
    std::vector<Index> indices;
    const size_t count = getFaceCount();
    for (Index i = 0; i < count; i++)
    {
      if (isBoundary(i))
        indices.push_back(i);
    }
    if (indices.size() == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Mesh has an empty boundary." << Alert::Raise;
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  FaceIterator MPIMesh::getInterface() const
  {
    std::vector<Index> indices;
    const size_t count = getFaceCount();
    for (Index i = 0; i < count; i++)
    {
      if (isInterface(i))
        indices.push_back(i);
    }
    if (indices.size() == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Mesh has an empty interface." << Alert::Raise;
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  bool MPIMesh::isInterface(Index faceIdx) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    bool local = false;
    const auto idx = getLocalIndex(getDimension() - 1, faceIdx);
    if (idx)
      local = shard.isInterface(*idx);
    return boost::mpi::all_reduce(comm, local, std::logical_or<bool>());
  }

  bool MPIMesh::isBoundary(Index faceIdx) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    bool local = false;
    const auto idx = getLocalIndex(getDimension() - 1, faceIdx);
    if (idx)
      local = shard.isBoundary(*idx);
    return boost::mpi::all_reduce(comm, local, std::logical_or<bool>());
  }

  Real MPIMesh::getMeasure(size_t d) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getMeasure(d);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getMeasure(size_t d, Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getMeasure(d, attr);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getMeasure(size_t d, const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getMeasure(d, attrs);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getVolume() const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getVolume();
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getVolume(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getVolume(attr);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getVolume(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getVolume(attrs);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getArea() const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getArea();
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getArea(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getArea(attr);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getArea(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getArea(attrs);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getPerimeter() const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getPerimeter();
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getPerimeter(Attribute attr) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getPerimeter(attr);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  Real MPIMesh::getPerimeter(const FlatSet<Attribute>& attrs) const
  {
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    Real local = shard.getPerimeter(attrs);
    return boost::mpi::all_reduce(comm, local, std::plus<Real>());
  }

  const PolytopeTransformation& MPIMesh::getPolytopeTransformation(size_t dimension, Index globalIdx) const
  {
    auto idx = getLocalIndex(dimension, globalIdx);
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    PolytopeTransformation* local;
    if (idx)
    {
      if (!shard.isGhost(dimension, *idx))
        local = shard.getPolytopeTransformation(dimension, *idx).copy();
    }
    auto res = boost::mpi::all_reduce(comm, local, [](auto const& a, auto const& b) { return a ? a : b; });
    assert(res);
    m_transformationIndex[dimension].write(
        [&](auto& obj)
        {
          assert(res);
          obj[globalIdx] = res;
        });
    return *res;
  }

  Polytope::Type MPIMesh::getGeometry(size_t dimension, Index globalIdx) const
  {
    auto idx = getLocalIndex(dimension, globalIdx);
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    boost::optional<Polytope::Type> local;
    if (idx)
    {
      if (!shard.isGhost(dimension, *idx))
        local = shard.getGeometry(dimension, *idx);
    }
    auto res = boost::mpi::all_reduce(comm, local, [](auto const& a, auto const& b) { return a ? a : b; });
    assert(res);
    return *res;
  }

  Attribute MPIMesh::getAttribute(size_t dimension, Index globalIdx) const
  {
    auto idx = getLocalIndex(dimension, globalIdx);
    const auto& comm = m_context.getCommunicator();
    const auto& shard = getShard();
    boost::optional<Attribute> local;
    if (idx)
    {
      if (!shard.isGhost(dimension, *idx))
        local = shard.getAttribute(dimension, *idx);
    }
    auto res = boost::mpi::all_reduce(comm, local, [](auto const& a, auto const& b) { return a ? a : b; });
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

  MPIMesh& MPIMesh::load(
    std::function<boost::filesystem::path(size_t)> filename, IO::FileFormat fmt)
  {
    const auto& comm = m_context.getCommunicator();
    return load(filename(comm.rank()), fmt);
  }

  MPIMesh& MPIMesh::load(
    const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    auto& shard = getShard();
    shard.load(filename, fmt);
    return *this;
  }

  void MPIMesh::save(
    std::function<boost::filesystem::path(size_t)> filename, IO::FileFormat fmt)
  {
    const auto& comm = m_context.getCommunicator();
    save(filename(comm.rank()), fmt);
  }

  void MPIMesh::save(
    const boost::filesystem::path& filename, IO::FileFormat fmt, size_t precison) const
  {
    const auto& shard = getShard();
    shard.save(filename, fmt);
  }

  Eigen::Map<const Math::PointVector> MPIMesh::getVertexCoordinates(Index globalIdx) const
  {
    auto idx = getLocalIndex(0, globalIdx);
    const auto& shard = getShard();
    const auto& comm = m_context.getCommunicator();
    Math::PointVector local;
    if (idx)
    {
      if (!shard.isGhost(0, *idx))
        local = shard.getVertexCoordinates(*idx);
    }
    assert(local.size() >= 0);
    assert(static_cast<size_t>(local.size()) == getSpaceDimension());
    auto res = boost::mpi::all_reduce(comm, local, [](auto const& a, auto const& b) { return a.size() > 0 ? a : b; });
    const auto& coords = (m_vertices[globalIdx] = res);
    return { coords.data(), static_cast<Eigen::Index>(coords.size()) };
  }
}

