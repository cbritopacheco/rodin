/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert/MemberFunctionException.h"

#include "Rodin/Variational/P1.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/MEDIT.h"

#include "Mesh.h"
#include "SubMesh.h"

#include "Polytope.h"
#include "PolytopeIterator.h"
#include "IsoparametricTransformation.h"

namespace Rodin::Geometry
{
  // ---- MeshBase ----------------------------------------------------------
  bool MeshBase::isSurface() const
  {
    return (getSpaceDimension() - 1 == getDimension());
  }

  // ---- Mesh<Context::Sequential> ----------------------------------------------
  Mesh<Context::Sequential>::Mesh(const Mesh& other)
    : m_sdim(other.m_sdim),
      m_vertices(other.m_vertices),
      m_connectivity(other.m_connectivity),
      m_attributeIndex(other.m_attributeIndex),
      m_attributes(other.m_attributes)
  {}

  Mesh<Context::Sequential>::Mesh(Mesh&& other)
    : m_sdim(std::move(other.m_sdim)),
      m_vertices(std::move(other.m_vertices)),
      m_connectivity(std::move(other.m_connectivity)),
      m_attributeIndex(std::move(other.m_attributeIndex)),
      m_transformationIndex(std::move(other.m_transformationIndex)),
      m_attributes(std::move(other.m_attributes))
  {}

  Mesh<Context::Sequential>& Mesh<Context::Sequential>::operator=(Mesh&& other)
  {
    Parent::operator=(std::move(other));
    m_sdim = std::move(other.m_sdim);
    m_vertices = std::move(other.m_vertices);
    m_connectivity = std::move(other.m_connectivity);
    m_attributeIndex = std::move(other.m_attributeIndex);
    m_transformationIndex = std::move(other.m_transformationIndex);
    m_attributes = std::move(other.m_attributes);
    return *this;
  }

  Mesh<Context::Sequential>::~Mesh()
  {
    for (auto& mt : m_transformationIndex)
    {
      mt.write(
          [](auto& obj)
          {
            for (PolytopeTransformation* ptr : obj)
              delete ptr;
          });
    }
  }

  Mesh<Context::Sequential>&
  Mesh<Context::Sequential>::load(const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    std::ifstream input(filename.c_str());
    if (!input)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to open " << filename << " for reading."
        << Alert::Raise;
    }
    switch (fmt)
    {
      case IO::FileFormat::MFEM:
      {
        IO::MeshLoader<IO::FileFormat::MFEM, Context::Sequential> loader(*this);
        loader.load(input);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::MeshLoader<IO::FileFormat::MEDIT, Context::Sequential> loader(*this);
        loader.load(input);
        break;
      }
      default:
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Loading from \"" << fmt << "\" format unsupported."
          << Alert::Raise;
        break;
      }
    }
    return *this;
  }

  void Mesh<Context::Sequential>::save(
      const boost::filesystem::path& filename,
      IO::FileFormat fmt, size_t precision) const
  {
    std::ofstream ofs(filename.c_str());
    if (!ofs)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to open " << filename << " for writing."
        << Alert::Raise;
    }
    ofs.precision(precision);
    switch (fmt)
    {
      case IO::FileFormat::MFEM:
      {
        IO::MeshPrinter<IO::FileFormat::MFEM, Context::Sequential> printer(*this);
        printer.print(ofs);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::MeshPrinter<IO::FileFormat::MEDIT, Context::Sequential> printer(*this);
        printer.print(ofs);
        break;
      }
      default:
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Saving to \"" << fmt << "\" format unsupported."
          << Alert::Raise;
      }
    }
    ofs.close();
  }

  SubMesh<Context::Sequential> Mesh<Context::Sequential>::keep(Attribute attr) const
  {
    return keep(FlatSet<Attribute>{attr});
  }

  SubMesh<Context::Sequential> Mesh<Context::Sequential>::keep(const FlatSet<Attribute>& attrs) const
  {
    const size_t D = getDimension();
    SubMesh<Context::Sequential>::Builder build;
    build.initialize(*this);
    for (Index i = 0; i < getCellCount(); i++)
    {
      if (attrs.count(getAttribute(D, i)))
      {
        build.include(D, i);
        for (size_t d = 1; d <= D - 1; d++)
        {
          const auto& inc = getConnectivity().getIncidence(D, d);
          if (inc.size() > 0)
            build.include(d, inc.at(i));
        }
      }
    }
    return build.finalize();
  }

  SubMesh<Context::Sequential> Mesh<Context::Sequential>::skin() const
  {
    const size_t D = getDimension();
    RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(D - 1, D);
    SubMesh<Context::Sequential>::Builder build;
    build.initialize(*this);
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      const Index i = it->getIndex();
      build.include(D - 1, i);
      for (size_t d = 1; d <= D - 2; d++)
      {
        const auto& inc = getConnectivity().getIncidence(D - 1, d);
        if (inc.size() > 0)
          build.include(d, inc.at(i));
      }
    }
    return build.finalize();
  }

  SubMesh<Context::Sequential> Mesh<Context::Sequential>::trim(Attribute attr) const
  {
    return trim(FlatSet<Attribute>{attr});
  }

  SubMesh<Context::Sequential> Mesh<Context::Sequential>::trim(const FlatSet<Attribute>& attrs) const
  {
    const size_t D = getDimension();
    SubMesh<Context::Sequential>::Builder build;
    build.initialize(*this);
    for (Index i = 0; i < getCellCount(); i++)
    {
      if (!attrs.count(getAttribute(D, i)))
      {
        build.include(D, i);
        for (size_t d = 1; d <= D - 1; d++)
        {
          const auto& inc = getConnectivity().getIncidence(D, d);
          if (inc.size() > 0)
            build.include(d, inc.at(i));
        }
      }
    }
    return build.finalize();
  }

  Mesh<Context::Sequential>& Mesh<Context::Sequential>::trace(
      const Map<std::pair<Attribute, Attribute>, Attribute>& interface, const FlatSet<Attribute>& attrs)
  {
    const size_t D = getDimension();
    RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(D - 1, D);
    const auto& conn = getConnectivity();
    for (auto it = getFace(); it; ++it)
    {
      if (attrs.size() == 0 || attrs.count(it->getAttribute()))
      {
        assert(it->getDimension() == D - 1);
        const auto& inc = conn.getIncidence({ D - 1, D }, it->getIndex());
        assert(inc.size() == 2);
        auto el1 = getCell(*inc.begin());
        auto el2 = getCell(*std::next(inc.begin()));
        auto find = interface.find({ el1->getAttribute(), el2->getAttribute() });
        if (find != interface.end())
        {
          setAttribute({ D - 1, it->getIndex() }, find->second);
        }
        else
        {
          find = interface.find({ el2->getAttribute(), el1->getAttribute() });
          if (find != interface.end())
            setAttribute({ D - 1, it->getIndex() }, find->second);
        }
      }
    }
    return *this;
  }

  Mesh<Context::Sequential>& Mesh<Context::Sequential>::scale(Scalar c)
  {
    m_vertices *= c;
    flush();
    return *this;
  }

#ifdef RODIN_USE_MPI
  Mesh<Context::MPI>
  Mesh<Context::Sequential>::parallelize(boost::mpi::communicator comm)
  {
    return Mesh<Context::MPI>(comm, *this);
  }
#endif

  Mesh<Context::Sequential>&
  Mesh<Context::Sequential>::setVertexCoordinates(Index idx, const Math::SpatialVector& coords)
  {
    m_vertices.col(idx) = coords;
    return *this;
  }

  Mesh<Context::Sequential>&
  Mesh<Context::Sequential>::setVertexCoordinates(Index idx, Scalar xi, size_t i)
  {
    m_vertices.col(idx).coeffRef(i) = xi;
    return *this;
  }

  Eigen::Map<const Math::SpatialVector> Mesh<Context::Sequential>::getVertexCoordinates(Index idx) const
  {
    const auto size = static_cast<Eigen::Index>(getSpaceDimension());
    return { getVertices().data() + getSpaceDimension() * idx, size };
  }

  const FlatSet<Attribute>& Mesh<Context::Sequential>::getAttributes(size_t d) const
  {
    return m_attributes[d];
  }

  size_t Mesh<Context::Sequential>::getDimension() const
  {
    return m_connectivity.getMeshDimension();
  }

  size_t Mesh<Context::Sequential>::getSpaceDimension() const
  {
    return m_sdim;
  }

  Mesh<Context::Sequential>& Mesh<Context::Sequential>::setPolytopeTransformation(
      const std::pair<size_t, Index> p, PolytopeTransformation* trans)
  {
    m_transformationIndex[p.first].write([&](auto& obj) { obj[p.second] = trans; });
    return *this;
  }

  PolytopeTransformation*
  Mesh<Context::Sequential>::getDefaultPolytopeTransformation(size_t dimension, Index idx) const
  {
    if (dimension == 0)
    {
      Variational::ScalarP1Element fe(Polytope::Type::Point);
      const size_t sdim = getSpaceDimension();
      Math::PointMatrix pm(sdim, 1);
      pm.col(0) = getVertexCoordinates(idx);
      return new IsoparametricTransformation(std::move(pm), std::move(fe));
    }
    else
    {
      auto g = getGeometry(dimension, idx);
      const size_t sdim = getSpaceDimension();
      const size_t n = Polytope::getVertexCount(g);
      Math::PointMatrix pm(sdim, n);
      const auto& polytope = getConnectivity().getPolytope(dimension, idx);
      assert(n == static_cast<size_t>(polytope.size()));
      for (const auto& v : polytope | boost::adaptors::indexed())
      {
        assert(sdim == static_cast<size_t>(getVertexCoordinates(v.value()).size()));
        pm.col(v.index()) = getVertexCoordinates(v.value());
      }
      Variational::ScalarP1Element fe(g);
      return new IsoparametricTransformation(std::move(pm), std::move(fe));
    }
  }

  const PolytopeTransformation&
  Mesh<Context::Sequential>::getPolytopeTransformation(size_t dimension, Index idx) const
  {
    assert(dimension < m_transformationIndex.size());
    if (m_transformationIndex[dimension].read().size() == 0)
    {
      m_transformationIndex[dimension].write(
          [&](auto& obj) { obj.resize(getPolytopeCount(dimension), nullptr); });
    }
    assert(0 < m_transformationIndex[dimension].read().size());
    assert(idx < m_transformationIndex[dimension].read().size());
    const auto& transPtr = m_transformationIndex[dimension].read()[idx];
    if (transPtr)
    {
      return *transPtr;
    }
    else
    {
      PolytopeTransformation* trans = getDefaultPolytopeTransformation(dimension, idx);
      m_transformationIndex[dimension].write(
          [&](auto& obj) { obj[idx] = trans; });
      return *trans;
    }
  }

  Scalar MeshBase::getVolume()
  {
    Scalar totalVolume = 0;
    for (auto it = getCell(); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Scalar MeshBase::getVolume(Attribute attr)
  {
    Scalar totalVolume = 0;
    for (auto it = getCell(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Scalar MeshBase::getPerimeter()
  {
    Scalar totalVolume = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Scalar MeshBase::getPerimeter(Attribute attr)
  {
    Scalar totalVolume = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  CCL MeshBase::ccl(
      std::function<bool(const Polytope&, const Polytope&)> p,
      size_t d,
      const FlatSet<Attribute>& attrs) const
  {
    FlatSet<Index> visited;
    visited.reserve(getPolytopeCount(d));
    std::deque<Index> searchQueue;
    std::deque<FlatSet<Index>> res;

    // Perform the labelling
    for (auto it = getPolytope(d); it; ++it)
    {
      const Index i = it->getIndex();
      if (!visited.count(i))
      {
        if (attrs.size() == 0 || attrs.count(it->getAttribute()))
        {
          res.push_back({});
          searchQueue.push_back(i);
        }
        while (searchQueue.size() > 0)
        {
          const Index idx = searchQueue.back();
          auto el = getPolytope(d, idx);
          searchQueue.pop_back();
          auto result = visited.insert(idx);
          const Boolean inserted = result.second;
          if (inserted)
          {
            res.back().insert(idx);
            for (auto adj = el->getAdjacent(); adj; ++adj)
            {
              if (p(*el, *adj))
              {
                if (attrs.size() == 0 || attrs.count(adj->getAttribute()))
                  searchQueue.push_back(adj->getIndex());
              }
            }
          }
        }
      }
    }
    return res;
  }

  size_t Mesh<Context::Sequential>::getPolytopeCount(size_t dimension) const
  {
    return m_connectivity.getCount(dimension);
  }

  size_t Mesh<Context::Sequential>::getPolytopeCount(Polytope::Type g) const
  {
    return m_connectivity.getCount(g);
  }

  FaceIterator Mesh<Context::Sequential>::getBoundary() const
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

  FaceIterator Mesh<Context::Sequential>::getInterface() const
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

  CellIterator Mesh<Context::Sequential>::getCell(Index idx) const
  {
    return CellIterator(*this, BoundedIndexGenerator(idx, getCellCount()));
  }

  FaceIterator Mesh<Context::Sequential>::getFace(Index idx) const
  {
    return FaceIterator(*this, BoundedIndexGenerator(idx, getFaceCount()));
  }

  VertexIterator Mesh<Context::Sequential>::getVertex(Index idx) const
  {
    return VertexIterator(*this, BoundedIndexGenerator(idx, getVertexCount()));
  }

  PolytopeIterator Mesh<Context::Sequential>::getPolytope(size_t dimension, Index idx) const
  {
    return PolytopeIterator(dimension, *this, BoundedIndexGenerator(idx, getPolytopeCount(dimension)));
  }

  bool Mesh<Context::Sequential>::isInterface(Index faceIdx) const
  {
    const size_t D = getDimension();
    RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(D - 1, D);
    const auto& conn = getConnectivity();
    assert(conn.getIncidence(D - 1, D).size());
    const auto& incidence = conn.getIncidence({D - 1, D}, faceIdx);
    assert(incidence.size() > 0);
    return incidence.size() > 1;
  }

  bool Mesh<Context::Sequential>::isBoundary(Index faceIdx) const
  {
    const size_t D = getDimension();
    RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(D - 1, D);
    const auto& conn = getConnectivity();
    assert(conn.getIncidence(D - 1, D).size());
    const auto& incidence = conn.getIncidence({D - 1, D}, faceIdx);
    assert(incidence.size() > 0);
    return incidence.size() == 1;
  }

  Polytope::Type Mesh<Context::Sequential>::getGeometry(size_t dimension, Index idx) const
  {
    return m_connectivity.getGeometry(dimension, idx);
  }

  Attribute Mesh<Context::Sequential>::getAttribute(size_t dimension, Index index) const
  {
    auto it = m_attributeIndex.find(dimension, index);
    if (it == m_attributeIndex.end(dimension))
      return RODIN_DEFAULT_POLYTOPE_ATTRIBUTE;
    else
      return it->second;
  }

  Mesh<Context::Sequential>&
  Mesh<Context::Sequential>::setAttribute(const std::pair<size_t, Index>& p, Attribute attr)
  {
    const auto [dimension, index] = p;
    m_attributeIndex.track(p, attr);
    m_attributes.at(dimension).insert(attr);
    return *this;
  }

  SubMeshBase& Mesh<Context::Sequential>::asSubMesh()
  {
    assert(isSubMesh());
    RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH();
    return static_cast<SubMesh<Context::Sequential>&>(*this);
  }

  const SubMeshBase& Mesh<Context::Sequential>::asSubMesh() const
  {
    assert(isSubMesh());
    RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH();
    return static_cast<const SubMesh<Context::Sequential>&>(*this);
  }

  Mesh<Context::Sequential> Mesh<Context::Sequential>::UniformGrid(Polytope::Type g, size_t h, size_t w)
  {
    Builder build;
    const size_t dim = Polytope::getGeometryDimension(g);
    build.initialize(dim).nodes(h * w);
    switch (g)
    {
      case Polytope::Type::Point:
      {
        return build.nodes(1).vertex({0}).finalize();
      }
      case Polytope::Type::Triangle:
      {
        assert(h * w >= 4);
        for (size_t i = 0; i < h; i++)
        {
          for (size_t j = 0; j < w; j++)
            build.vertex({ static_cast<Scalar>(j), static_cast<Scalar>(i) });
        }

        build.reserve(dim, 2 * (h - 1) * (w - 1));
        for (size_t i = 0; i < h - 1; i++)
        {
          for (size_t j = 0; j < w - 1; j++)
          {
            build.polytope(g, { i * w + j, i * w + j + 1, (i + 1) * w + j })
                 .polytope(g, { i * w + j + 1, (i + 1) * w + j + 1, (i + 1) * w + j });
          }
        }
        return build.finalize();
      }
      default:
      {
        assert(false);
        return build.nodes(0).finalize();
      }
    };
  }
}

