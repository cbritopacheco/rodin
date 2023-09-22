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

  // ---- Mesh<Context::Serial> ----------------------------------------------
  Mesh<Context::Serial>&
  Mesh<Context::Serial>::load(const boost::filesystem::path& filename, IO::FileFormat fmt)
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
        IO::MeshLoader<IO::FileFormat::MFEM, Context::Serial> loader(*this);
        loader.load(input);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::MeshLoader<IO::FileFormat::MEDIT, Context::Serial> loader(*this);
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

  void Mesh<Context::Serial>::save(
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
        IO::MeshPrinter<IO::FileFormat::MFEM, Context::Serial> printer(*this);
        printer.print(ofs);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::MeshPrinter<IO::FileFormat::MEDIT, Context::Serial> printer(*this);
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

  SubMesh<Context::Serial> Mesh<Context::Serial>::keep(Attribute attr) const
  {
    return keep(FlatSet<Attribute>{attr});
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::keep(const FlatSet<Attribute>& attrs) const
  {
    const size_t D = getDimension();
    SubMesh<Context::Serial>::Builder build;
    build.initialize(*this);
    for (Index i = 0; i < getElementCount(); i++)
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

  SubMesh<Context::Serial> Mesh<Context::Serial>::skin() const
  {
    const size_t D = getDimension();
    SubMesh<Context::Serial>::Builder build;
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

  SubMesh<Context::Serial> Mesh<Context::Serial>::trim(Attribute attr) const
  {
    return trim(FlatSet<Attribute>{attr});
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::trim(const FlatSet<Attribute>& attrs) const
  {
    const size_t D = getDimension();
    SubMesh<Context::Serial>::Builder build;
    build.initialize(*this);
    for (Index i = 0; i < getElementCount(); i++)
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

  Mesh<Context::Serial>& Mesh<Context::Serial>::trace(
      const Map<std::pair<Attribute, Attribute>, Attribute>& tmap)
  {
    const size_t D = getDimension();
    for (auto it = getFace(); it; ++it)
    {
      assert(it->getDimension() == D - 1);
      const auto& inc = getConnectivity().getIncidence({ D - 1, D }, it->getIndex());
      assert(inc.size() == 2);
    }
    return *this;
  }

  Mesh<Context::Serial>& Mesh<Context::Serial>::scale(Scalar c)
  {
    m_vertices *= c;
    flush();
    return *this;
  }

#ifdef RODIN_USE_MPI
  Mesh<Context::MPI>
  Mesh<Context::Serial>::parallelize(boost::mpi::communicator comm)
  {
    return Mesh<Context::MPI>(comm, *this);
  }
#endif

  Eigen::Map<const Math::SpatialVector> Mesh<Context::Serial>::getVertexCoordinates(Index idx) const
  {
    const auto size = static_cast<Eigen::Index>(getSpaceDimension());
    return { getVertices().data() + getSpaceDimension() * idx, size };
  }

  const FlatSet<Attribute>& Mesh<Context::Serial>::getAttributes(size_t d) const
  {
    return m_attributes[d];
  }

  size_t Mesh<Context::Serial>::getDimension() const
  {
    return m_connectivity.getMeshDimension();
  }

  size_t Mesh<Context::Serial>::getSpaceDimension() const
  {
    return m_sdim;
  }

  const PolytopeTransformation&
  Mesh<Context::Serial>::getPolytopeTransformation(size_t dimension, Index idx) const
  {
    auto it = m_transformationIndex.find(dimension, idx);
    if (it != m_transformationIndex.end(dimension))
    {
      assert(m_transformationIndex.at(dimension, idx));
      return *it->second;
    }
    else
    {
      if (dimension == 0)
      {
        Variational::ScalarP1Element fe(Polytope::Type::Point);
        const size_t sdim = getSpaceDimension();
        Math::SpatialMatrix pm(sdim, 1);
        pm.col(0) = getVertexCoordinates(idx);
        auto trans =
          std::unique_ptr<PolytopeTransformation>(
              new IsoparametricTransformation(pm, std::move(fe)));
        auto p = m_transformationIndex.insert(it, { dimension, idx }, std::move(trans));
        return *p->second;
      }
      else
      {
        auto g = getGeometry(dimension, idx);
        const size_t sdim = getSpaceDimension();
        const size_t n = Polytope::getVertexCount(g);
        Math::SpatialMatrix pm(sdim, n);
        const auto& polytope = getConnectivity().getPolytope(dimension, idx);
        assert(n == static_cast<size_t>(polytope.size()));
        for (const auto& v : polytope | boost::adaptors::indexed())
        {
          assert(sdim == static_cast<size_t>(getVertexCoordinates(v.value()).size()));
          pm.col(v.index()) = getVertexCoordinates(v.value());
        }
        Variational::ScalarP1Element fe(g);
        auto trans =
          std::unique_ptr<PolytopeTransformation>(
              new IsoparametricTransformation(std::move(pm), std::move(fe)));
        auto p = m_transformationIndex.insert(it, { dimension, idx }, std::move(trans));
        return *p->second;
      }
    }
  }

  Scalar MeshBase::getVolume()
  {
    Scalar totalVolume = 0;
    for (auto it = getElement(); !it.end(); ++it)
      totalVolume += it->getMeasure();
    return totalVolume;
  }

  Scalar MeshBase::getVolume(Attribute attr)
  {
    Scalar totalVolume = 0;
    for (auto it = getElement(); !it.end(); ++it)
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

  // std::deque<std::set<int>> MeshBase::ccl(
  //     std::function<bool(const Element&, const Element&)> p) const
  // {
  //   std::set<int> visited;
  //   std::deque<int> searchQueue;
  //   std::deque<std::set<int>> res;

  //   // Perform the labelling
  //   assert(false);
  //   // for (int i = 0; i < count<Element>(); i++)
  //   // {
  //   //   if (!visited.count(i))
  //   //   {
  //   //     res.push_back({});
  //   //     searchQueue.push_back(i);
  //   //     while (searchQueue.size() > 0)
  //   //     {
  //   //       int el = searchQueue.back();
  //   //       searchQueue.pop_back();
  //   //       auto result = visited.insert(el);
  //   //       bool inserted = result.second;
  //   //       if (inserted)
  //   //       {
  //   //         res.back().insert(el);
  //   //         for (int n : get<Element>(el).adjacent())
  //   //         {
  //   //           if (p(get<Element>(el), get<Element>(n)))
  //   //           {
  //   //             searchQueue.push_back(n);
  //   //           }
  //   //         }
  //   //       }
  //   //     }
  //   //   }
  //   // }
  //   return res;
  // }

  size_t Mesh<Context::Serial>::getCount(size_t dimension) const
  {
    return m_connectivity.getCount(dimension);
  }

  size_t Mesh<Context::Serial>::getCount(Polytope::Type g) const
  {
    return m_connectivity.getCount(g);
  }

  FaceIterator Mesh<Context::Serial>::getBoundary() const
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
        << "Mesh has empty boundary." << Alert::Raise;
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  FaceIterator Mesh<Context::Serial>::getInterface() const
  {
    std::vector<Index> indices;
    const size_t count = getFaceCount();
    for (Index i = 0; i < count; i++)
    {
      if (isInterface(i))
        indices.push_back(i);
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  ElementIterator Mesh<Context::Serial>::getElement(Index idx) const
  {
    return ElementIterator(*this, BoundedIndexGenerator(idx, getElementCount()));
  }

  FaceIterator Mesh<Context::Serial>::getFace(Index idx) const
  {
    return FaceIterator(*this, BoundedIndexGenerator(idx, getFaceCount()));
  }

  VertexIterator Mesh<Context::Serial>::getVertex(Index idx) const
  {
    return VertexIterator(*this, BoundedIndexGenerator(idx, getVertexCount()));
  }

  PolytopeIterator Mesh<Context::Serial>::getPolytope(size_t dimension, Index idx) const
  {
    return PolytopeIterator(dimension, *this, BoundedIndexGenerator(idx, getCount(dimension)));
  }

  bool Mesh<Context::Serial>::isInterface(Index faceIdx) const
  {
    const size_t D = getDimension();
    const auto& incidence = getConnectivity().getIncidence({D - 1, D}, faceIdx);
    assert(incidence.size() > 0);
    return incidence.size() > 1;
  }

  bool Mesh<Context::Serial>::isBoundary(Index faceIdx) const
  {
    const size_t D = getDimension();
    const auto& conn = getConnectivity();
    if (conn.getIncidence(D - 1, D).size() == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << Alert::Notation::Incidence(D - 1, D)
        << " has not been computed and is required to use this function."
        << Alert::Raise;
    }
    assert(conn.getIncidence(D - 1, D).size());
    const auto& incidence = conn.getIncidence({D - 1, D}, faceIdx);
    assert(incidence.size() > 0);
    return incidence.size() == 1;
  }

  Polytope::Type Mesh<Context::Serial>::getGeometry(size_t dimension, Index idx) const
  {
    return m_connectivity.getGeometry(dimension, idx);
  }

  Attribute Mesh<Context::Serial>::getAttribute(size_t dimension, Index index) const
  {
    auto it = m_attributeIndex.find(dimension, index);
    if (it == m_attributeIndex.end(dimension))
      return RODIN_DEFAULT_POLYTOPE_ATTRIBUTE;
    else
      return it->second;
  }

  Mesh<Context::Serial>&
  Mesh<Context::Serial>::setAttribute(const std::pair<size_t, Index>& p, Attribute attr)
  {
    const auto [dimension, index] = p;
    m_attributeIndex.track(p, attr);
    m_attributes.at(dimension).insert(attr);
    return *this;
  }

  SubMeshBase& Mesh<Context::Serial>::asSubMesh()
  {
    assert(isSubMesh());
    if (!isSubMesh())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "This instance of Mesh is not a SubMesh (isSubMesh() == false). "
        << "Downcasting to SubMesh is ill-defined."
        << Alert::Raise;
    }
    return static_cast<SubMesh<Context::Serial>&>(*this);
  }

  const SubMeshBase& Mesh<Context::Serial>::asSubMesh() const
  {
    assert(isSubMesh());
    if (!isSubMesh())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "This instance of Mesh is not a SubMesh (isSubMesh() == false). "
        << "Downcasting to SubMesh is ill-defined."
        << Alert::Raise;
    }
    return static_cast<const SubMesh<Context::Serial>&>(*this);
  }

  Mesh<Context::Serial> Mesh<Context::Serial>::UniformGrid(Polytope::Type g, size_t h, size_t w)
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

