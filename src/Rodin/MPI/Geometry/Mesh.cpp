#include <boost/dynamic_bitset.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeTransformation.h"

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
}

