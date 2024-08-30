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

#include "ForwardDecls.h"

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
        IO::MeshLoader<IO::FileFormat::MFEM, Context> loader(*this);
        loader.load(input);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::MeshLoader<IO::FileFormat::MEDIT, Context> loader(*this);
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
        IO::MeshPrinter<IO::FileFormat::MFEM, Context> printer(*this);
        printer.print(ofs);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::MeshPrinter<IO::FileFormat::MEDIT, Context> printer(*this);
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
    SubMesh<Context>::Builder build;
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
    SubMesh<Context>::Builder build;
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
    SubMesh<Context>::Builder build;
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

  Mesh<Context::Sequential>& Mesh<Context::Sequential>::scale(Real c)
  {
    m_vertices *= c;
    flush();
    return *this;
  }

  Mesh<Context::Sequential>&
  Mesh<Context::Sequential>::setVertexCoordinates(Index idx, const Math::SpatialVector<Real>& coords)
  {
    m_vertices.col(idx) = coords;
    return *this;
  }

  Mesh<Context::Sequential>&
  Mesh<Context::Sequential>::setVertexCoordinates(Index idx, Real xi, size_t i)
  {
    m_vertices.col(idx).coeffRef(i) = xi;
    return *this;
  }

  Eigen::Map<const Math::SpatialVector<Real>> Mesh<Context::Sequential>::getVertexCoordinates(Index idx) const
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
      Variational::RealP1Element fe(Polytope::Type::Point);
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
      Variational::RealP1Element fe(g);
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

  Real MeshBase::getVolume() const
  {
    Real totalVolume = 0;
    for (auto it = getCell(); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Real MeshBase::getVolume(Attribute attr) const
  {
    Real totalVolume = 0;
    for (auto it = getCell(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real MeshBase::getVolume(const FlatSet<Attribute>& attrs) const
  {
    Real totalVolume = 0;
    for (auto it = getCell(); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real MeshBase::getPerimeter() const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
      totalPerimeter += it->getMeasure();
    return totalPerimeter;
  }

  Real MeshBase::getPerimeter(Attribute attr) const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalPerimeter += it->getMeasure();
    }
    return totalPerimeter;
  }

  Real MeshBase::getPerimeter(const FlatSet<Attribute>& attrs) const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalPerimeter += it->getMeasure();
    }
    return totalPerimeter;
  }

  Real MeshBase::getArea() const
  {
    Real totalArea = 0;
    for (auto it = getFace(); !it.end(); ++it)
      totalArea += it->getMeasure();
    return totalArea;
  }

  Real MeshBase::getArea(Attribute attr) const
  {
    Real totalArea = 0;
    for (auto it = getFace(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalArea += it->getMeasure();
    }
    return totalArea;
  }

  Real MeshBase::getArea(const FlatSet<Attribute>& attrs) const
  {
    Real totalArea = 0;
    for (auto it = getFace(); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalArea += it->getMeasure();
    }
    return totalArea;
  }

  Real MeshBase::getMeasure(size_t d) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Real MeshBase::getMeasure(size_t d, Attribute attr) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real MeshBase::getMeasure(size_t d, const FlatSet<Attribute>& attrs) const
  {
    Real totalMeasure = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalMeasure += it->getMeasure();
    }
    return totalMeasure;
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
    return static_cast<SubMesh<Context>&>(*this);
  }

  const SubMeshBase& Mesh<Context::Sequential>::asSubMesh() const
  {
    assert(isSubMesh());
    RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH();
    return static_cast<const SubMesh<Context>&>(*this);
  }

  Mesh<Context::Sequential> Mesh<Context::Sequential>::UniformGrid(Polytope::Type g, const Array<size_t>& dimensions)
  {
    Builder build;
    const size_t dim = Polytope::getGeometryDimension(g);
    switch (g)
    {
      case Polytope::Type::Point:
      {
        return build.nodes(1).vertex({0}).finalize();
      }
      case Polytope::Type::Triangle:
      {
        assert(dimensions.size() == 2);
        const size_t w = dimensions.coeff(0);
        const size_t h = dimensions.coeff(1);
        build.initialize(dim).nodes(w * h);
        assert(w * h >= 4);
        for (size_t j = 0; j < h; j++)
        {
          for (size_t i = 0; i < w; i++)
            build.vertex({ static_cast<Real>(i), static_cast<Real>(j) });
        }

        build.reserve(dim, 2 * (h - 1) * (w - 1));
        for (size_t i = 0; i < w - 1; i++)
        {
          for (size_t j = 0; j < h - 1; j++)
          {
            build.polytope(g, { i + j * w, (i + 1) + j * w , i + (j + 1) * w })
                 .polytope(g, { (i + 1) + j * w, (i + 1) + (j + 1) * w, i + (j + 1) * w });
          }
        }
        return build.finalize();
      }
      case Polytope::Type::Tetrahedron:
      {
        assert(dimensions.size() == 3);
        const size_t width = dimensions.coeff(0);
        const size_t height = dimensions.coeff(1);
        const size_t depth = dimensions.coeff(2);
        assert(width * height * depth >= 8);
        build.initialize(dim)
             .nodes(width * height * depth + (width - 1) * (height - 1) * (depth - 1))
             .reserve(dim, 10 * (width - 1) * (height - 1) * (depth - 1));

        for (size_t k = 0; k < depth; ++k)
        {
          for (size_t j = 0; j < height; ++j)
          {
            for (size_t i = 0; i < width; ++i)
            {
              build.vertex({
                  static_cast<Real>(i),
                  static_cast<Real>(j),
                  static_cast<Real>(k) });
            }
          }
        }

        for (size_t k = 0; k < depth - 1; ++k)
        {
          for (size_t j = 0; j < height - 1; ++j)
          {
            for (size_t i = 0; i < width - 1; ++i)
            {
              build.vertex({
                  static_cast<Real>(i + 0.5),
                  static_cast<Real>(j + 0.5),
                  static_cast<Real>(k + 0.5) });
            }
          }
        }

        for (size_t i = 0; i < width - 1; ++i)
        {
          for (size_t j = 0; j < height - 1; ++j)
          {
            for (size_t k = 0; k < depth - 1; ++k)
            {
              const Index c =
                  i + (width - 1) * j + (width - 1) * (height - 1) * k
                    + (width - 1) + width * (height - 1) + width * height * (depth - 1) + 1;
              build.polytope(g, // Front-left
                       { i + width * j + width * height * k,
                        (i + 1) + width * j + width * height * k,
                         i + width * (j + 1) + width * height * k,
                         i + width * j + width * height * (k + 1) })
                   .polytope(g, // Front-right
                       { (i + 1) + width * j + width * height * (k + 1),
                          i + width * j + width * height * (k + 1),
                          c,
                          (i + 1) + width * j + width * height * k })
                   .polytope(g, // Left-top
                       { i + width * (j + 1) + width * height * (k + 1),
                         c,
                         i + width * j + width * height * (k + 1),
                         i + width * (j + 1) + width * height * k })
                   .polytope(g, // Top-left
                       { i + width * j + width * height * (k + 1),
                        (i + 1) + width * j + width * height * (k + 1),
                         i + width * (j + 1) + width * height * (k + 1),
                         c })
                   .polytope(g, // Right-bottom
                       { c,
                        (i + 1) + width * j + width * height * k,
                         (i + 1) + width * (j + 1) + width * height * k,
                         (i + 1) + width * j + width * height * (k + 1) })
                   .polytope(g, // Bottom-right
                       { (i + 1) + width * j + width * height * k,
                         i + width * (j + 1) + width * height * k,
                         (i + 1) + width * (j + 1) + width * height * k,
                          c
                         })
                   .polytope(g, // Back-left
                       {  i + width * (j + 1) + width * height * k,
                          (i + 1) + width * (j + 1) + width * height * k,
                          c,
                          i + width * (j + 1) + width * height * (k + 1) })
                   .polytope(g, // Back-right
                        { (i + 1) + width * (j + 1) + width * height * (k + 1),
                           i + width * (j + 1) + width * height * (k + 1),
                           (i + 1) + width * j + width * height * (k + 1),
                          (i + 1) + width * (j + 1) + width * height * k })
                   .polytope(g, // Front fill
                       { (i + 1) + width * j + width * height * k,
                          i + width * (j + 1) + width * height * k,
                          c,
                          i + width * j + width * height * (k + 1) })
                   .polytope(g, // Back fill
                       { (i + 1) + width * j + width * height * (k + 1),
                          (i + 1) + width * (j + 1) + width * height * k,
                          c,
                          i + width * (j + 1) + width * height * (k + 1) })
                   ;
            }
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

  std::optional<Point> Mesh<Context::Sequential>::inclusion(const Point& p) const
  {
    const auto& polytope = p.getPolytope();
    if (!polytope.getMesh().isSubMesh())
    {
      return {};
    }
    const auto& submesh = polytope.getMesh().asSubMesh();
    const auto& ancestors = submesh.getAncestors();
    const size_t d = polytope.getDimension();
    Index i = polytope.getIndex();
    i = submesh.getPolytopeMap(d).left.at(i);
    auto it = ancestors.begin();
    while (it != ancestors.end())
    {
      if (it->get() == *this)
      {
        auto pit = this->getPolytope(d, i);
        std::unique_ptr<Polytope> parentPolytope(pit.release());
        return Point(
            std::move(*parentPolytope),
            this->getPolytopeTransformation(d, i),
            std::cref(p.getReferenceCoordinates()),
            p.getPhysicalCoordinates());
      }
      else if (it->get().isSubMesh())
      {
        const auto& parentMesh = it->get().asSubMesh();
        i = parentMesh.getPolytopeMap(d).left.at(i);
      }
      else
      {
        // Invalid inclusion.
        // The SubMesh where the Point belongs to is not a descendant of this Mesh.
        return {};
      }
      ++it;
    }
    // Invalid inclusion.
    // The SubMesh where the Point belongs to is not a descendant of this Mesh.
    return {};
  }
}

