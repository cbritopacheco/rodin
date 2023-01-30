/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"
#include "Rodin/IO/MeshLoader.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Mesh.h"
#include "SubMesh.h"

#include "Simplex.h"
#include "SimplexIterator.h"


namespace Rodin::Geometry
{
  // ---- MeshBase ----------------------------------------------------------
  bool MeshBase::isSurface() const
  {
    return (getSpaceDimension() - 1 == getDimension());
  }

  std::set<Attribute> MeshBase::getAttributes() const
  {
    return std::set<Attribute>(
        getHandle().attributes.begin(), getHandle().attributes.end());
  }

  std::set<Attribute> MeshBase::getBoundaryAttributes() const
  {
    return std::set<Attribute>(
        getHandle().bdr_attributes.begin(), getHandle().bdr_attributes.end());
  }

  // ---- Mesh<Context::Serial> ----------------------------------------------
  size_t Mesh<Context::Serial>::getDimension() const
  {
    return m_dim;
  }

  size_t Mesh<Context::Serial>::getSpaceDimension() const
  {
    return m_sdim;
  }

  Mesh<Context::Serial>& Mesh<Context::Serial>::scale(double c)
  {
    mfem::Vector vs;
    getHandle().GetVertices(vs);
    vs *= c;
    getHandle().SetVertices(vs);
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
      case IO::FileFormat::GMSH:
      {
        IO::MeshPrinter<IO::FileFormat::GMSH, Context::Serial> printer(*this);
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

  MeshBase& MeshBase::displace(const Variational::GridFunctionBase& u)
  {
    assert(u.getFiniteElementSpace().getVectorDimension() == getSpaceDimension());
    getHandle().MoveNodes(u.getHandle());
    return *this;
  }

  double
  MeshBase::getMaximumDisplacement(const Variational::GridFunctionBase& u)
  {
    double res;
    getHandle().CheckDisplacements(u.getHandle(), res);
    return res;
  }

  double MeshBase::getVolume()
  {
    double totalVolume = 0;
    for (auto it = getElement(); !it.end(); ++it)
      totalVolume += it->getVolume();
    return totalVolume;
  }

  double MeshBase::getVolume(Attribute attr)
  {
    double totalVolume = 0;
    for (auto it = getElement(); !it.end(); ++it)
    {
      if (it->getAttribute() == attr)
        totalVolume += it->getVolume();
    }
    return totalVolume;
  }

  double MeshBase::getPerimeter()
  {
    double totalVolume = 0;
    for (auto it = getBoundary(); !it.end(); ++it)
      totalVolume += it->getVolume();
    return totalVolume;
  }

  double MeshBase::getPerimeter(Attribute attr)
  {
    double totalVolume = 0;
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

  Mesh<Context::Serial>::Mesh(mfem::Mesh&& mesh)
    : m_dim(mesh.Dimension()), m_sdim(mesh.SpaceDimension()), m_impl(std::move(mesh))
  {
    for (int i = 0; i < getHandle().GetNBE(); i++)
      m_f2b[getHandle().GetBdrElementEdgeIndex(i)] = i;
  }

  size_t Mesh<Context::Serial>::getCount(size_t dimension) const
  {
    if (dimension == getDimension())
    {
      return m_impl.GetNE();
    }
    else if (dimension == getDimension() - 1)
    {
      return m_impl.GetNumFaces();
    }
    else if (dimension == 0)
    {
      return m_impl.GetNV();
    }
    else
    {
      assert(false);
    }
    return 0;
  }

  FaceIterator Mesh<Context::Serial>::getBoundary() const
  {
    std::vector<Index> indices;
    indices.reserve(getHandle().GetNBE());
    for (int i = 0; i < getHandle().GetNBE(); i++)
    {
      int idx = getHandle().GetBdrFace(i);
      if (!getHandle().FaceIsInterior(idx))
        indices.push_back(idx);
    }
    return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
  }

  FaceIterator Mesh<Context::Serial>::getInterface() const
  {
    std::vector<Index> indices;
    indices.reserve(getHandle().GetNumFaces());
    for (int idx = 0; idx < getHandle().GetNumFaces(); idx++)
    {
      if (getHandle().FaceIsInterior(idx))
        indices.push_back(idx);
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

  SimplexIterator Mesh<Context::Serial>::getSimplex(size_t dimension, Index idx) const
  {
    return SimplexIterator(dimension, *this, BoundedIndexGenerator(idx, getCount(dimension)));
  }

  bool Mesh<Context::Serial>::isInterface(Index faceIdx) const
  {
    return getHandle().FaceIsInterior(faceIdx);
  }

  bool Mesh<Context::Serial>::isBoundary(Index faceIdx) const
  {
    return !getHandle().FaceIsInterior(faceIdx);
  }

  Attribute Mesh<Context::Serial>::getAttribute(size_t dimension, Index index) const
  {
    if (dimension == getDimension())
    {
      return getHandle().GetAttribute(index);
    }
    else if (dimension == getDimension() - 1)
    {
      auto it = m_f2b.find(index);
      if (it != m_f2b.end())
      {
        return getHandle().GetBdrAttribute(it->second);
      }
      else
      {
        return RODIN_DEFAULT_SIMPLEX_ATTRIBUTE;
      }
    }
    else if (dimension == 0)
    {
      return RODIN_DEFAULT_SIMPLEX_ATTRIBUTE;
    }
    else
    {
      assert(false);
      return RODIN_DEFAULT_SIMPLEX_ATTRIBUTE;
    }
  }

  Mesh<Context::Serial>& Mesh<Context::Serial>
  ::setAttribute(size_t dimension, Index index, Attribute attr)
  {
    if (dimension == getDimension())
    {
      getHandle().SetAttribute(index, attr);
    }
    else if (dimension == getDimension() - 1)
    {
      auto it = m_f2b.find(index);
      if (it != m_f2b.end())
      {
        getHandle().SetBdrAttribute(it->second, attr);
      }
      else
      {
        assert(false);
      }
    }
    else
    {
    }
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
      case IO::FileFormat::GMSH:
      {
        IO::MeshLoader<IO::FileFormat::GMSH, Context::Serial> loader(*this);
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
        Alert::Exception()
          << "Loading from \"" << fmt << "\" format unsupported."
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
    SubMesh<Context::Serial> res(*this);
    std::set<Index> indices;
    for (Index i = 0; i < getCount(getDimension()); i++)
    {
      if (attrs.count(getAttribute(getDimension(), i)))
        indices.insert(i);
    }
    res.initialize(getDimension(), getSpaceDimension())
       .include(getDimension(), indices)
       .finalize();
    return res;
  }

  Mesh<Context::Serial>::Builder
  Mesh<Context::Serial>::initialize(size_t dim, size_t sdim)
  {
    m_dim = dim;
    m_sdim = sdim;
    Mesh<Context::Serial>::Builder build;
    build.setReference(*this);
    return build;
  }

  SubMesh<Context::Serial> Mesh<Context::Serial>::skin() const
  {
    assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!
    SubMesh<Context::Serial> res(*this);
    std::set<Index> indices;
    for (auto it = getBoundary(); !it.end(); ++it)
      indices.insert(it->getIndex());
    res.initialize(getDimension() - 1, getSpaceDimension())
       .include(getDimension() - 1, indices)
       .finalize();
    return res;
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

  mfem::Mesh& Mesh<Context::Serial>::getHandle()
  {
    return m_impl;
  }

  const mfem::Mesh& Mesh<Context::Serial>::getHandle() const
  {
    return m_impl;
  }
}

