/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/MEDIT.h"

#include "Mesh.h"
#include "SubMesh.h"

#include "Simplex.h"
#include "SimplexIterator.h"
#include "IsoparametricTransformation.h"

namespace Rodin::Geometry
{
  // ---- MeshBase ----------------------------------------------------------
  bool MeshBase::isSurface() const
  {
    return (getSpaceDimension() - 1 == getDimension());
  }

  // ---- Mesh<Context::Serial> ----------------------------------------------
  Eigen::Map<const Math::Vector> Mesh<Context::Serial>::getVertexCoordinates(Index idx) const
  {
    const auto size = static_cast<Eigen::Index>(getSpaceDimension());
    return { getVertices().data() + getSpaceDimension() * idx, size };
  }

  const std::set<Attribute>& Mesh<Context::Serial>::getAttributes() const
  {
    return m_attributes[getDimension()];
  }

  const std::set<Attribute>& Mesh<Context::Serial>::getBoundaryAttributes() const
  {
    return m_attributes[getDimension()];
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
    if (m_transformations.isTracked(dimension, idx))
    {
      assert(m_transformations.at(dimension, idx));
      return *m_transformations.at(dimension, idx);
    }
    else
    {
      auto g = getGeometry(dimension, idx);
      const size_t sdim = getSpaceDimension();
      const size_t n = Polytope::getVertexCount(g);
      Math::Matrix pm(sdim, n);
      const auto& polytope = getConnectivity().getPolytope(dimension, idx);
      assert(n == static_cast<size_t>(polytope.size()));
      for (const auto& v : polytope | boost::adaptors::indexed())
      {
        assert(sdim == static_cast<size_t>(getVertexCoordinates(v.value()).size()));
        pm.col(v.index()) = getVertexCoordinates(v.value());
      }
      Variational::P1Element fe(g);
      auto trans =
        std::unique_ptr<PolytopeTransformation>(
            new IsoparametricTransformation(std::move(pm), std::move(fe)));
      m_transformations.track(dimension, idx, std::move(trans));
      return *m_transformations.at(dimension, idx);
    }
  }

  Mesh<Context::Serial>& Mesh<Context::Serial>::scale(Scalar c)
  {
    m_vertices *= c;
    return *this;
  }

  void Mesh<Context::Serial>::save(
      const boost::filesystem::path& filename,
      IO::FileFormat fmt, size_t precision) const
  {
    std::ofstream ofs(filename.c_str());
    if (!ofs)
    {
      Alert::Exception()
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
        Alert::Exception()
          << "Saving to \"" << fmt << "\" format unsupported."
          << Alert::Raise;
      }
    }
  }

  Scalar MeshBase::getVolume()
  {
    Scalar totalVolume = 0;
    for (auto it = getElement(); !it.end(); ++it)
      totalVolume += it->getVolume();
    return totalVolume;
  }

  Scalar MeshBase::getVolume(Attribute attr)
  {
    Scalar totalVolume = 0;
    for (auto it = getElement(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getVolume();
    }
    return totalVolume;
  }

  Scalar MeshBase::getPerimeter()
  {
    Scalar totalVolume = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
      totalVolume += it->getVolume();
    return totalVolume;
  }

  Scalar MeshBase::getPerimeter(Attribute attr)
  {
    Scalar totalVolume = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getVolume();
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

  // ---- Mesh<Serial> ------------------------------------------------------
#ifdef RODIN_USE_MPI
  Mesh<Context::MPI>
  Mesh<Context::Serial>::parallelize(boost::mpi::communicator comm)
  {
    return Mesh<Context::MPI>(comm, *this);
  }
#endif

  size_t Mesh<Context::Serial>::getCount(size_t dimension) const
  {
    return m_connectivity.getCount(dimension);
  }

  size_t Mesh<Context::Serial>::getCount(Polytope::Geometry g) const
  {
    return m_connectivity.getCount(g);
  }

  FaceIterator Mesh<Context::Serial>::getBoundary() const
  {
    assert(false);
    // std::vector<Index> indices;
    // indices.reserve(getHandle().GetNBE());
    // for (int i = 0; i < getHandle().GetNBE(); i++)
    // {
    //   int idx = getHandle().GetBdrFace(i);
    //   if (!getHandle().FaceIsInterior(idx))
    //     indices.push_back(idx);
    // }
    // return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  FaceIterator Mesh<Context::Serial>::getInterface() const
  {
    assert(false);
    // std::vector<Index> indices;
    // indices.reserve(getHandle().GetNumFaces());
    // for (int idx = 0; idx < getHandle().GetNumFaces(); idx++)
    // {
    //   if (getHandle().FaceIsInterior(idx))
    //     indices.push_back(idx);
    // }
    // return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
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
    assert(false);
    // return getHandle().FaceIsInterior(faceIdx);
  }

  bool Mesh<Context::Serial>::isBoundary(Index faceIdx) const
  {
    assert(false);
    // return !getHandle().FaceIsInterior(faceIdx);
  }

  Polytope::Geometry Mesh<Context::Serial>::getGeometry(size_t dimension, Index idx) const
  {
    return m_connectivity.getGeometry(dimension, idx);
  }

  Attribute Mesh<Context::Serial>::getAttribute(size_t dimension, Index index) const
  {
    auto it = m_attrs.find(dimension, index);
    if (it == m_attrs.end(dimension))
      return RODIN_DEFAULT_POLYTOPE_ATTRIBUTE;
    else
      return it->second;
  }

  Mesh<Context::Serial>&
  Mesh<Context::Serial>::setAttribute(size_t dimension, Index index, Attribute attr)
  {
    m_attrs.track(dimension, index, attr);
    return *this;
  }

  Mesh<Context::Serial>&
  Mesh<Context::Serial>::load(const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    mfem::named_ifgzstream input(filename.c_str());
    if (!input)
    {
      Alert::Exception()
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
        Alert::Exception() << "Loading from \"" << fmt << "\" format unsupported."
                           << Alert::Raise;
        break;
      }
    }

    return *this;
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::keep(Attribute attr)
  {
    return keep(std::set<Attribute>{attr});
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::keep(const std::set<Attribute>& attrs)
  {
    assert(false);
    // SubMesh<Context::Serial> res(*this);
    // std::set<Index> indices;
    // for (Index i = 0; i < getCount(getDimension()); i++)
    // {
    //   if (attrs.count(getAttribute(getDimension(), i)))
    //     indices.insert(i);
    // }
    // res.initialize(getSpaceDimension()).include(indices).finalize();
    // return res;
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::skin() const
  {
    assert(false);
    // SubMesh<Context::Serial> res(*this);
    // std::set<Index> indices;
    // for (auto it = getBoundary(); !it.end(); ++it)
    //   indices.insert(it->getIndex());
    // res.initialize(getSpaceDimension()).include(indices).finalize();
    // return res;
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::trim(Attribute attr)
  {
    return trim(std::set<Attribute>{attr});
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::trim(const std::set<Attribute>& attrs)
  {
    std::set<Attribute> complement = getAttributes();
    for (const auto& a : attrs)
      complement.erase(a);
    return keep(complement);
  }
}

