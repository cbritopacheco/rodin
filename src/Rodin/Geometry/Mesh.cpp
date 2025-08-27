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
#include "Rodin/IO/EnSight6.h"

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

  // ---- Mesh<Context::Local> ----------------------------------------------
  Mesh<Context::Local>::Mesh(const Mesh& other)
    : m_sdim(other.m_sdim),
      m_vertices(other.m_vertices),
      m_connectivity(other.m_connectivity),
      m_attributeIndex(other.m_attributeIndex),
      m_attributes(other.m_attributes)
  {}

  Mesh<Context::Local>::Mesh(Mesh&& other)
    : m_sdim(std::move(other.m_sdim)),
      m_vertices(std::move(other.m_vertices)),
      m_connectivity(std::move(other.m_connectivity)),
      m_attributeIndex(std::move(other.m_attributeIndex)),
      m_transformationIndex(std::move(other.m_transformationIndex)),
      m_attributes(std::move(other.m_attributes))
  {}

  Mesh<Context::Local>& Mesh<Context::Local>::operator=(Mesh&& other)
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

  Mesh<Context::Local>::~Mesh()
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

  Mesh<Context::Local>&
  Mesh<Context::Local>::load(const boost::filesystem::path& filename, IO::FileFormat fmt)
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

  void Mesh<Context::Local>::save(
      const boost::filesystem::path& filename, IO::FileFormat fmt) const
  {
    std::ofstream ofs(filename.c_str());
    if (!ofs)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to open " << filename << " for writing."
        << Alert::Raise;
    }
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
      case IO::FileFormat::ENSIGHT6:
      {
        IO::MeshPrinter<IO::FileFormat::ENSIGHT6, Context> printer(*this);
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

  SubMesh<Context::Local> Mesh<Context::Local>::keep(Attribute attr) const
  {
    return keep(FlatSet<Attribute>{attr});
  }

  SubMesh<Context::Local> Mesh<Context::Local>::keep(const FlatSet<Attribute>& attrs) const
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

  SubMesh<Context::Local> Mesh<Context::Local>::skin() const
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

  SubMesh<Context::Local> Mesh<Context::Local>::trim(Attribute attr) const
  {
    return trim(FlatSet<Attribute>{attr});
  }

  SubMesh<Context::Local> Mesh<Context::Local>::trim(const FlatSet<Attribute>& attrs) const
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


  Mesh<Context::Local>& Mesh<Context::Local>::trace(
      const Map<Attribute, Attribute>& interface, const FlatSet<Attribute>& attrs)
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
        if (inc.size() == 1)
        {
          auto el = getCell(*inc.begin());
          auto find = interface.find(el->getAttribute());
          if (find != interface.end())
            setAttribute({ D - 1, it->getIndex() }, find->second);
        }
      }
    }
    return *this;
  }

  Mesh<Context::Local>& Mesh<Context::Local>::trace(
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
        if (inc.size() == 2)
        {
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
    }
    return *this;
  }

  Mesh<Context::Local>& Mesh<Context::Local>::scale(Real c)
  {
    m_vertices *= c;
    flush();
    return *this;
  }

  Mesh<Context::Local>&
  Mesh<Context::Local>::setVertexCoordinates(Index idx, const Math::SpatialVector<Real>& coords)
  {
    m_vertices.col(idx) = coords;
    return *this;
  }

  Mesh<Context::Local>&
  Mesh<Context::Local>::setVertexCoordinates(Index idx, Real xi, size_t i)
  {
    m_vertices.col(idx).coeffRef(i) = xi;
    return *this;
  }

  Eigen::Map<const Math::SpatialVector<Real>> Mesh<Context::Local>::getVertexCoordinates(Index idx) const
  {
    const auto size = static_cast<Eigen::Index>(getSpaceDimension());
    return { m_vertices.data() + getSpaceDimension() * idx, size };
  }

  const FlatSet<Attribute>& Mesh<Context::Local>::getAttributes(size_t d) const
  {
    return m_attributes[d];
  }

  size_t Mesh<Context::Local>::getDimension() const
  {
    return m_connectivity.getMeshDimension();
  }

  size_t Mesh<Context::Local>::getSpaceDimension() const
  {
    return m_sdim;
  }

  Mesh<Context::Local>& Mesh<Context::Local>::setPolytopeTransformation(
      const std::pair<size_t, Index> p, PolytopeTransformation* trans)
  {
    m_transformationIndex[p.first].write([&](auto& obj) { obj[p.second] = trans; });
    return *this;
  }

  PolytopeTransformation*
  Mesh<Context::Local>::getDefaultPolytopeTransformation(size_t dimension, Index idx) const
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
      const size_t n = Polytope::Traits(g).getVertexCount();
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
  Mesh<Context::Local>::getPolytopeTransformation(size_t dimension, Index idx) const
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

  Real Mesh<Context::Local>::getVolume() const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(3); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Real Mesh<Context::Local>::getVolume(Attribute attr) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(3); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real Mesh<Context::Local>::getVolume(const FlatSet<Attribute>& attrs) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(3); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real Mesh<Context::Local>::getPerimeter() const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
      totalPerimeter += it->getMeasure();
    return totalPerimeter;
  }

  Real Mesh<Context::Local>::getPerimeter(Attribute attr) const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalPerimeter += it->getMeasure();
    }
    return totalPerimeter;
  }

  Real Mesh<Context::Local>::getPerimeter(const FlatSet<Attribute>& attrs) const
  {
    Real totalPerimeter = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalPerimeter += it->getMeasure();
    }
    return totalPerimeter;
  }

  Real Mesh<Context::Local>::getArea() const
  {
    Real totalArea = 0;
    for (auto it = getPolytope(2); !it.end(); ++it)
      totalArea += it->getMeasure();
    return totalArea;
  }

  Real Mesh<Context::Local>::getArea(Attribute attr) const
  {
    Real totalArea = 0;
    for (auto it = getPolytope(2); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalArea += it->getMeasure();
    }
    return totalArea;
  }

  Real Mesh<Context::Local>::getArea(const FlatSet<Attribute>& attrs) const
  {
    Real totalArea = 0;
    for (auto it = getPolytope(2); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalArea += it->getMeasure();
    }
    return totalArea;
  }

  Real Mesh<Context::Local>::getMeasure(size_t d) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Real Mesh<Context::Local>::getMeasure(size_t d, Attribute attr) const
  {
    Real totalVolume = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getMeasure();
    }
    return totalVolume;
  }

  Real Mesh<Context::Local>::getMeasure(size_t d, const FlatSet<Attribute>& attrs) const
  {
    Real totalMeasure = 0;
    for (auto it = getPolytope(d); !it.end(); ++it)
    {
      if (attrs.contains(it->getAttribute()))
        totalMeasure += it->getMeasure();
    }
    return totalMeasure;
  }

  size_t Mesh<Context::Local>::getPolytopeCount(size_t dimension) const
  {
    return m_connectivity.getCount(dimension);
  }

  size_t Mesh<Context::Local>::getPolytopeCount(Polytope::Type g) const
  {
    return m_connectivity.getCount(g);
  }

  FaceIterator Mesh<Context::Local>::getBoundary() const
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

  FaceIterator Mesh<Context::Local>::getInterface() const
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

  CellIterator Mesh<Context::Local>::getCell() const
  {
    return getCell(0);
  }

  FaceIterator Mesh<Context::Local>::getFace() const
  {
    return getFace(0);
  }

  VertexIterator Mesh<Context::Local>::getVertex() const
  {
    return getVertex(0);
  }

  PolytopeIterator Mesh<Context::Local>::getPolytope(size_t dimension) const
  {
    return getPolytope(dimension, 0);
  }

  CellIterator Mesh<Context::Local>::getCell(Index idx) const
  {
    return CellIterator(*this, BoundedIndexGenerator(idx, getCellCount()));
  }

  FaceIterator Mesh<Context::Local>::getFace(Index idx) const
  {
    return FaceIterator(*this, BoundedIndexGenerator(idx, getFaceCount()));
  }

  VertexIterator Mesh<Context::Local>::getVertex(Index idx) const
  {
    return VertexIterator(*this, BoundedIndexGenerator(idx, getVertexCount()));
  }

  PolytopeIterator Mesh<Context::Local>::getPolytope(size_t dimension, Index idx) const
  {
    return PolytopeIterator(dimension, *this, BoundedIndexGenerator(idx, getPolytopeCount(dimension)));
  }

  bool Mesh<Context::Local>::isInterface(Index faceIdx) const
  {
    const size_t D = getDimension();
    RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(D - 1, D);
    const auto& conn = getConnectivity();
    assert(conn.getIncidence(D - 1, D).size());
    const auto& incidence = conn.getIncidence({D - 1, D}, faceIdx);
    assert(incidence.size() > 0);
    return incidence.size() > 1;
  }

  bool Mesh<Context::Local>::isBoundary(Index faceIdx) const
  {
    const size_t D = getDimension();
    RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(D - 1, D);
    const auto& conn = getConnectivity();
    assert(conn.getIncidence(D - 1, D).size());
    const auto& incidence = conn.getIncidence({D - 1, D}, faceIdx);
    assert(incidence.size() > 0);
    return incidence.size() == 1;
  }

  Polytope::Type Mesh<Context::Local>::getGeometry(size_t dimension, Index idx) const
  {
    return m_connectivity.getGeometry(dimension, idx);
  }

  Attribute Mesh<Context::Local>::getAttribute(size_t dimension, Index index) const
  {
    auto it = m_attributeIndex.find(dimension, index);
    if (it == m_attributeIndex.end(dimension))
      return RODIN_DEFAULT_POLYTOPE_ATTRIBUTE;
    else
      return it->second;
  }

  Mesh<Context::Local>&
  Mesh<Context::Local>::setAttribute(const std::pair<size_t, Index>& p, Attribute attr)
  {
    const auto [dimension, index] = p;
    m_attributeIndex.track(p, attr);
    m_attributes.at(dimension).insert(attr);
    return *this;
  }

  SubMeshBase& Mesh<Context::Local>::asSubMesh()
  {
    assert(isSubMesh());
    RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH();
    return static_cast<SubMesh<Context>&>(*this);
  }

  const SubMeshBase& Mesh<Context::Local>::asSubMesh() const
  {
    assert(isSubMesh());
    RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH();
    return static_cast<const SubMesh<Context>&>(*this);
  }

  Optional<Point> MeshBase::inclusion(const Point& p) const
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
        return Point(*this->getPolytope(d, i), p.getReferenceCoordinates(), p.getPhysicalCoordinates());
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

  Mesh<Context::Local> Mesh<Context::Local>::UniformGrid(Polytope::Type g, const Array<size_t>& dimensions)
  {
    Builder build;
    const size_t dim = Polytope::Traits(g).getDimension();
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
      case Polytope::Type::Quadrilateral:
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

        build.reserve(dim, (h - 1) * (w - 1));
        for (size_t i = 0; i < w - 1; i++)
        {
          for (size_t j = 0; j < h - 1; j++)
          {
            build.polytope(g, {
                i + j * w, (i + 1) + j * w , i + (j + 1) * w,  (i + 1) + (j + 1) * w });
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
      case Polytope::Type::Wedge:
      {
        assert(dimensions.size() == 3);
        const size_t w = dimensions.coeff(0);
        const size_t h = dimensions.coeff(1);
        const size_t d = dimensions.coeff(2);
        assert(w >= 2 && h >= 2 && d >= 2);
        build.initialize(dim).nodes(w * h * d);
        for (size_t k = 0; k < d; k++)
        {
          for (size_t j = 0; j < h; j++)
          {
            for (size_t i = 0; i < w; i++)
            {
              build.vertex({ static_cast<Real>(i),
                             static_cast<Real>(j),
                             static_cast<Real>(k) });
            }
          }
        }
        build.reserve(dim, 2 * (w - 1) * (h - 1) * (d - 1));
        for (size_t k = 0; k < d - 1; k++)
        {
          for (size_t j = 0; j < h - 1; j++)
          {
            for (size_t i = 0; i < w - 1; i++)
            {
              const Index v0 = i + j * w + k * (w * h);
              const Index v1 = (i + 1) + j * w + k * (w * h);
              const Index v2 = i + (j + 1) * w + k * (w * h);
              const Index v3 = (i + 1) + (j + 1) * w + k * (w * h);
              const Index v0p = v0 + w * h;
              const Index v1p = v1 + w * h;
              const Index v2p = v2 + w * h;
              const Index v3p = v3 + w * h;
              build.polytope(g, { v0, v1, v2, v0p, v1p, v2p });
              build.polytope(g, { v1, v3, v2, v1p, v3p, v2p });
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

  Mesh<Context::Local>
  Mesh<Context::Local>::Box(Polytope::Type faceType, const Array<size_t>& n)
  {
    const size_t dim  = Polytope::Traits(faceType).getDimension();
    const size_t sdim = dim + 1;
    assert(static_cast<size_t>(n.size()) == sdim);

    Builder build;
    build.initialize(sdim);

    switch (faceType)
    {
      // 0D in 1D: endpoints only
      case Polytope::Type::Point:
      {
        const size_t Nx = n.coeff(0);
        assert(Nx >= 1);
        build.nodes(2);
        build.vertex({ Real(0) });
        build.vertex({ Real(Nx) });
        return build.finalize();
      }

      // 1D in 2D: rectangle boundary only
      case Polytope::Type::Segment:
      {
        const size_t Nx = n.coeff(0), Ny = n.coeff(1);
        assert(Nx >= 1 && Ny >= 1);

        build.nodes(2 * Nx + 2 * Ny);

        for (size_t i = 0; i <= Nx; ++i) { build.vertex({ Real(i), Real(0) }); }
        for (size_t j = 1; j <= Ny; ++j) { build.vertex({ Real(Nx), Real(j) }); }
        for (size_t i = Nx; i-- > 0; )  { build.vertex({ Real(i),  Real(Ny) }); }
        for (size_t j = Ny; j-- > 1; )  { build.vertex({ Real(0),  Real(j)  }); }

        build.reserve(1, 2 * Nx + 2 * Ny);

        // bottom
        for (size_t i = 0; i < Nx; ++i)
        {
          build.polytope(Polytope::Type::Segment, { Index(i), Index(i + 1) });
        }
        // right
        if (Ny >= 1)
        {
          build.polytope(Polytope::Type::Segment, { Index(Nx), Index(Nx + 1) });
          for (size_t j = 0; j + 1 < Ny; ++j)
          {
            build.polytope(Polytope::Type::Segment, { Index(Nx + 1 + j), Index(Nx + 1 + j + 1) });
          }
        }
        // top
        if (Nx >= 1)
        {
          build.polytope(Polytope::Type::Segment, { Index(Nx + (Ny ? Ny : 0)), Index(Nx + Ny + 1) });
          for (size_t t = 0; t + 1 < Nx; ++t)
          {
            build.polytope(Polytope::Type::Segment, { Index(Nx + Ny + 1 + t), Index(Nx + Ny + 1 + t + 1) });
          }
        }
        // left
        if (Ny >= 1)
        {
          if (Ny > 1)
          {
            build.polytope(Polytope::Type::Segment, { Index(Nx + Ny + (Nx ? Nx : 0)), Index(Nx + Ny + Nx + 1) });
            for (size_t l = 0; l + 1 < Ny - 1; ++l)
            {
              build.polytope(Polytope::Type::Segment, { Index(Nx + Ny + Nx + 1 + l), Index(Nx + Ny + Nx + 1 + l + 1) });
            }
            build.polytope(Polytope::Type::Segment, { Index(Nx + Ny + Nx + (Ny - 1)), Index(0) });
          }
          else
          {
            build.polytope(Polytope::Type::Segment, { Index(Nx + Ny + (Nx ? Nx : 0)), Index(0) });
          }
        }

        return build.finalize();
      }

      // 2D in 3D: triangle surface mesh
      case Polytope::Type::Triangle:
      {
        const size_t Nx = n.coeff(0), Ny = n.coeff(1), Nz = n.coeff(2);
        assert(Nx >= 1 && Ny >= 1 && Nz >= 1);

        const size_t Lx = Nx + 1, Ly = Ny + 1;

        // vertex counts by blocks
        const size_t cntZ0 = Lx * Ly;                                  // z=0
        const size_t cntZN = Lx * Ly;                                  // z=Nz
        const size_t cntX0 = (Nz >= 2 ? (Nz - 1) * Ly : 0);            // x=0, k=1..Nz-1
        const size_t cntXN = (Nz >= 2 ? (Nz - 1) * Ly : 0);            // x=Nx, k=1..Nz-1
        const size_t cntY0 = (Nx >= 2 && Nz >= 2 ? (Nx - 1) * (Nz - 1) : 0); // y=0,  i=1..Nx-1,k=1..Nz-1
        const size_t cntYN = cntY0;                                    // y=Ny

        // offsets
        const size_t offZ0 = 0;
        const size_t offZN = offZ0 + cntZ0;
        const size_t offX0 = offZN + cntZN;
        const size_t offXN = offX0 + cntX0;
        const size_t offY0 = offXN + cntXN;
        const size_t offYN = offY0 + cntY0;

        build.nodes(offYN + cntYN);

        // emit vertices: z=0
        for (size_t j = 0; j < Ly; ++j)
        {
          for (size_t i = 0; i < Lx; ++i)
          {
            build.vertex({ Real(i), Real(j), Real(0) });
          }
        }
        // z=Nz
        for (size_t j = 0; j < Ly; ++j)
        {
          for (size_t i = 0; i < Lx; ++i)
          {
            build.vertex({ Real(i), Real(j), Real(Nz) });
          }
        }
        // x=0, k=1..Nz-1
        if (Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t j = 0; j < Ly; ++j)
            {
              build.vertex({ Real(0), Real(j), Real(k) });
            }
          }
        }
        // x=Nx, k=1..Nz-1
        if (Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t j = 0; j < Ly; ++j)
            {
              build.vertex({ Real(Nx), Real(j), Real(k) });
            }
          }
        }
        // y=0, i=1..Nx-1, k=1..Nz-1
        if (Nx >= 2 && Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t i = 1; i <= Nx - 1; ++i)
            {
              build.vertex({ Real(i), Real(0), Real(k) });
            }
          }
        }
        // y=Ny, i=1..Nx-1, k=1..Nz-1
        if (Nx >= 2 && Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t i = 1; i <= Nx - 1; ++i)
            {
              build.vertex({ Real(i), Real(Ny), Real(k) });
            }
          }
        }

        build.reserve(2, 4 * (Ny * Nz + Nx * Nz + Nx * Ny));

        // z = 0  (-z outward) flip winding
        for (size_t j = 0; j < Ny; ++j)
        {
          for (size_t i = 0; i < Nx; ++i)
          {
            Index a = Index(offZ0 + (i + 0) + (j + 0) * Lx);
            Index b = Index(offZ0 + (i + 1) + (j + 0) * Lx);
            Index c = Index(offZ0 + (i + 1) + (j + 1) * Lx);
            Index d = Index(offZ0 + (i + 0) + (j + 1) * Lx);
            build.polytope(Polytope::Type::Triangle, { a, c, b });
            build.polytope(Polytope::Type::Triangle, { a, d, c });
          }
        }
        // z = Nz (+z outward)
        for (size_t j = 0; j < Ny; ++j)
        {
          for (size_t i = 0; i < Nx; ++i)
          {
            Index a = Index(offZN + (i + 0) + (j + 0) * Lx);
            Index b = Index(offZN + (i + 0) + (j + 1) * Lx);
            Index c = Index(offZN + (i + 1) + (j + 1) * Lx);
            Index d = Index(offZN + (i + 1) + (j + 0) * Lx);
            build.polytope(Polytope::Type::Triangle, { a, b, c });
            build.polytope(Polytope::Type::Triangle, { a, c, d });
          }
        }
        // x = 0  (-x outward)
        if (Nz >= 1)
        {
          for (size_t j = 0; j < Ny; ++j)
          {
            for (size_t k = 0; k < Nz; ++k)
            {
              Index a = (k == 0)      ? Index(offZ0 + 0 + (j + 0) * Lx) : Index(offX0 + (k - 1) * Ly + (j + 0));
              Index b = (k == 0)      ? Index(offZ0 + 0 + (j + 1) * Lx) : Index(offX0 + (k - 1) * Ly + (j + 1));
              Index c = (k + 1 == Nz) ? Index(offZN + 0 + (j + 1) * Lx) : Index(offX0 + (k + 0) * Ly + (j + 1));
              Index d = (k + 1 == Nz) ? Index(offZN + 0 + (j + 0) * Lx) : Index(offX0 + (k + 0) * Ly + (j + 0));
              build.polytope(Polytope::Type::Triangle, { a, b, c });
              build.polytope(Polytope::Type::Triangle, { a, c, d });
            }
          }
        }
        // x = Nx (+x outward)
        if (Nz >= 1)
        {
          for (size_t j = 0; j < Ny; ++j)
          {
            for (size_t k = 0; k < Nz; ++k)
            {
              Index a = (k == 0)      ? Index(offZ0 + Nx + (j + 0) * Lx) : Index(offXN + (k - 1) * Ly + (j + 0));
              Index b = (k + 1 == Nz) ? Index(offZN + Nx + (j + 0) * Lx) : Index(offXN + (k + 0) * Ly + (j + 0));
              Index c = (k + 1 == Nz) ? Index(offZN + Nx + (j + 1) * Lx) : Index(offXN + (k + 0) * Ly + (j + 1));
              Index d = (k == 0)      ? Index(offZ0 + Nx + (j + 1) * Lx) : Index(offXN + (k - 1) * Ly + (j + 1));
              build.polytope(Polytope::Type::Triangle, { a, b, c });
              build.polytope(Polytope::Type::Triangle, { a, c, d });
            }
          }
        }
        // y = 0  (-y outward)
        if (Nz >= 1)
        {
          for (size_t i = 0; i < Nx; ++i)
          {
            for (size_t k = 0; k < Nz; ++k)
            {
              Index a = (k == 0)      ? Index(offZ0 + (i + 0) + 0 * Lx)
                                      : (i == 0 ? Index(offX0 + (k - 1) * Ly + 0)
                                                : (i == Nx ? Index(offXN + (k - 1) * Ly + 0)
                                                           : Index(offY0 + (k - 1) * (Lx - 2) + (i - 1))));
              Index b = (k + 1 == Nz) ? Index(offZN + (i + 0) + 0 * Lx)
                                      : (i == 0 ? Index(offX0 + (k + 0) * Ly + 0)
                                                : (i == Nx ? Index(offXN + (k + 0) * Ly + 0)
                                                           : Index(offY0 + (k + 0) * (Lx - 2) + (i - 1))));
              Index c = (k + 1 == Nz) ? Index(offZN + (i + 1) + 0 * Lx)
                                      : (i + 1 == Nx ? Index(offXN + (k + 0) * Ly + 0)
                                                     : Index(offY0 + (k + 0) * (Lx - 2) + (i + 0)));
              Index d = (k == 0)      ? Index(offZ0 + (i + 1) + 0 * Lx)
                                      : (i + 1 == Nx ? Index(offXN + (k - 1) * Ly + 0)
                                                     : Index(offY0 + (k - 1) * (Lx - 2) + (i + 0)));
              build.polytope(Polytope::Type::Triangle, { a, b, c });
              build.polytope(Polytope::Type::Triangle, { a, c, d });
            }
          }
        }
        // y = Ny (+y outward)
        if (Nz >= 1)
        {
          for (size_t i = 0; i < Nx; ++i)
          {
            for (size_t k = 0; k < Nz; ++k)
            {
              Index a = (k == 0)      ? Index(offZ0 + (i + 0) + Ny * Lx)
                                      : (i == 0 ? Index(offX0 + (k - 1) * Ly + Ny)
                                                : (i == Nx ? Index(offXN + (k - 1) * Ly + Ny)
                                                           : Index(offYN + (k - 1) * (Lx - 2) + (i - 1))));
              Index b = (k == 0)      ? Index(offZ0 + (i + 1) + Ny * Lx)
                                      : (i + 1 == Nx ? Index(offXN + (k - 1) * Ly + Ny)
                                                     : Index(offYN + (k - 1) * (Lx - 2) + (i + 0)));
              Index c = (k + 1 == Nz) ? Index(offZN + (i + 1) + Ny * Lx)
                                      : (i + 1 == Nx ? Index(offXN + (k + 0) * Ly + Ny)
                                                     : Index(offYN + (k + 0) * (Lx - 2) + (i + 0)));
              Index d = (k + 1 == Nz) ? Index(offZN + (i + 0) + Ny * Lx)
                                      : (i == 0 ? Index(offX0 + (k + 0) * Ly + Ny)
                                                : (i == Nx ? Index(offXN + (k + 0) * Ly + Ny)
                                                           : Index(offYN + (k + 0) * (Lx - 2) + (i - 1))));
              build.polytope(Polytope::Type::Triangle, { a, b, c });
              build.polytope(Polytope::Type::Triangle, { a, c, d });
            }
          }
        }

        return build.finalize();
      }

      // 2D in 3D: quadrilateral surface mesh
      case Polytope::Type::Quadrilateral:
      {
        const size_t Nx = n.coeff(0), Ny = n.coeff(1), Nz = n.coeff(2);
        assert(Nx >= 1 && Ny >= 1 && Nz >= 1);

        const size_t Lx = Nx + 1, Ly = Ny + 1;

        // vertex counts
        const size_t cntZ0 = Lx * Ly;
        const size_t cntZN = Lx * Ly;
        const size_t cntX0 = (Nz >= 2 ? (Nz - 1) * Ly : 0);
        const size_t cntXN = (Nz >= 2 ? (Nz - 1) * Ly : 0);
        const size_t cntY0 = (Nx >= 2 && Nz >= 2 ? (Nx - 1) * (Nz - 1) : 0);
        const size_t cntYN = cntY0;

        // offsets
        const size_t offZ0 = 0;
        const size_t offZN = offZ0 + cntZ0;
        const size_t offX0 = offZN + cntZN;
        const size_t offXN = offX0 + cntX0;
        const size_t offY0 = offXN + cntXN;
        const size_t offYN = offY0 + cntY0;

        build.nodes(offYN + cntYN);

        // emit vertices
        for (size_t j = 0; j < Ly; ++j)
        {
          for (size_t i = 0; i < Lx; ++i)
          {
            build.vertex({ Real(i), Real(j), Real(0) });
          }
        }
        for (size_t j = 0; j < Ly; ++j)
        {
          for (size_t i = 0; i < Lx; ++i)
          {
            build.vertex({ Real(i), Real(j), Real(Nz) });
          }
        }
        if (Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t j = 0; j < Ly; ++j)
            {
              build.vertex({ Real(0), Real(j), Real(k) });
            }
          }
        }
        if (Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t j = 0; j < Ly; ++j)
            {
              build.vertex({ Real(Nx), Real(j), Real(k) });
            }
          }
        }
        if (Nx >= 2 && Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t i = 1; i <= Nx - 1; ++i)
            {
              build.vertex({ Real(i), Real(0), Real(k) });
            }
          }
        }
        if (Nx >= 2 && Nz >= 2)
        {
          for (size_t k = 1; k <= Nz - 1; ++k)
          {
            for (size_t i = 1; i <= Nx - 1; ++i)
            {
              build.vertex({ Real(i), Real(Ny), Real(k) });
            }
          }
        }

        build.reserve(2, 2 * (Ny * Nz + Nx * Nz + Nx * Ny));

        // z = 0  (−z outward) order (0,0)(1,0)(0,1)(1,1)
        for (size_t j = 0; j < Ny; ++j)
        {
          for (size_t i = 0; i < Nx; ++i)
          {
            build.polytope(Polytope::Type::Quadrilateral, {
              Index(offZ0 + (i + 0) + (j + 0) * Lx),
              Index(offZ0 + (i + 1) + (j + 0) * Lx),
              Index(offZ0 + (i + 0) + (j + 1) * Lx),
              Index(offZ0 + (i + 1) + (j + 1) * Lx)
            });
          }
        }

        // z = Nz (+z outward) order (0,0)(0,1)(1,0)(1,1)
        for (size_t j = 0; j < Ny; ++j)
        {
          for (size_t i = 0; i < Nx; ++i)
          {
            build.polytope(Polytope::Type::Quadrilateral, {
              Index(offZN + (i + 0) + (j + 0) * Lx),
              Index(offZN + (i + 0) + (j + 1) * Lx),
              Index(offZN + (i + 1) + (j + 0) * Lx),
              Index(offZN + (i + 1) + (j + 1) * Lx)
            });
          }
        }

        // x = 0  (−x outward) param (y,z): (0,0)(1,0)(0,1)(1,1)
        for (size_t j = 0; j < Ny; ++j)
        {
          for (size_t k = 0; k < Nz; ++k)
          {
            build.polytope(Polytope::Type::Quadrilateral, {
              Index((k == 0)      ? (offZ0 + 0 + (j + 0) * Lx) : (offX0 + (k - 1) * Ly + (j + 0))),
              Index((k == 0)      ? (offZ0 + 0 + (j + 1) * Lx) : (offX0 + (k - 1) * Ly + (j + 1))),
              Index((k + 1 == Nz) ? (offZN + 0 + (j + 0) * Lx) : (offX0 + (k + 0) * Ly + (j + 0))),
              Index((k + 1 == Nz) ? (offZN + 0 + (j + 1) * Lx) : (offX0 + (k + 0) * Ly + (j + 1)))
            });
          }
        }

        // x = Nx (+x outward) param (y,z): (0,0)(1,0)(0,1)(1,1)
        for (size_t j = 0; j < Ny; ++j)
        {
          for (size_t k = 0; k < Nz; ++k)
          {
            build.polytope(Polytope::Type::Quadrilateral, {
              Index((k == 0)      ? (offZ0 + Nx + (j + 0) * Lx) : (offXN + (k - 1) * Ly + (j + 0))),
              Index((k + 1 == Nz) ? (offZN + Nx + (j + 0) * Lx) : (offXN + (k + 0) * Ly + (j + 0))),
              Index((k == 0)      ? (offZ0 + Nx + (j + 1) * Lx) : (offXN + (k - 1) * Ly + (j + 1))),
              Index((k + 1 == Nz) ? (offZN + Nx + (j + 1) * Lx) : (offXN + (k + 0) * Ly + (j + 1)))
            });
          }
        }

        // y = 0  (−y outward) param (x,z): (0,0)(1,0)(0,1)(1,1)
        for (size_t i = 0; i < Nx; ++i)
        {
          for (size_t k = 0; k < Nz; ++k)
          {
            build.polytope(Polytope::Type::Quadrilateral, {
              Index((k == 0)      ? (offZ0 + (i + 0) + 0 * Lx)
                                   : (i == 0 ? (offX0 + (k - 1) * Ly + 0)
                                             : (i == Nx ? (offXN + (k - 1) * Ly + 0)
                                                        : (offY0 + (k - 1) * (Lx - 2) + (i - 1))))),
              Index((k == 0)      ? (offZ0 + (i + 1) + 0 * Lx)
                                   : (i + 1 == Nx ? (offXN + (k - 1) * Ly + 0)
                                                  : (offY0 + (k - 1) * (Lx - 2) + (i + 0)))),
              Index((k + 1 == Nz) ? (offZN + (i + 0) + 0 * Lx)
                                   : (i == 0 ? (offX0 + (k + 0) * Ly + 0)
                                             : (i == Nx ? (offXN + (k + 0) * Ly + 0)
                                                        : (offY0 + (k + 0) * (Lx - 2) + (i - 1))))),
              Index((k + 1 == Nz) ? (offZN + (i + 1) + 0 * Lx)
                                   : (i + 1 == Nx ? (offXN + (k + 0) * Ly + 0)
                                                  : (offY0 + (k + 0) * (Lx - 2) + (i + 0))))
            });
          }
        }

        // y = Ny (+y outward) param (x,z): (0,0)(1,0)(0,1)(1,1)
        for (size_t i = 0; i < Nx; ++i)
        {
          for (size_t k = 0; k < Nz; ++k)
          {
            build.polytope(Polytope::Type::Quadrilateral, {
              Index((k == 0)      ? (offZ0 + (i + 0) + Ny * Lx)
                                   : (i == 0 ? (offX0 + (k - 1) * Ly + Ny)
                                             : (i == Nx ? (offXN + (k - 1) * Ly + Ny)
                                                        : (offYN + (k - 1) * (Lx - 2) + (i - 1))))),
              Index((k == 0)      ? (offZ0 + (i + 1) + Ny * Lx)
                                   : (i + 1 == Nx ? (offXN + (k - 1) * Ly + Ny)
                                                  : (offYN + (k - 1) * (Lx - 2) + (i + 0)))),
              Index((k + 1 == Nz) ? (offZN + (i + 0) + Ny * Lx)
                                   : (i == 0 ? (offX0 + (k + 0) * Ly + Ny)
                                             : (i == Nx ? (offXN + (k + 0) * Ly + Ny)
                                                        : (offYN + (k + 0) * (Lx - 2) + (i - 1))))),
              Index((k + 1 == Nz) ? (offZN + (i + 1) + Ny * Lx)
                                   : (i + 1 == Nx ? (offXN + (k + 0) * Ly + Ny)
                                                  : (offYN + (k + 0) * (Lx - 2) + (i + 0))))
            });
          }
        }

        return build.finalize();
      }

      default:
      {
        assert(false && "Unsupported face type for Box");
        return build.nodes(0).finalize();
      }
    }
  }
}

