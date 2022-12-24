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

#include "Element.h"
#include "SimplexIterator.h"


namespace Rodin::Geometry
{
   // ---- MeshBase ----------------------------------------------------------
   size_t MeshBase::getSpaceDimension() const
   {
      return getHandle().SpaceDimension();
   }

   size_t MeshBase::getDimension() const
   {
      return getHandle().Dimension();
   }

   bool MeshBase::isSurface() const
   {
      return (getSpaceDimension() - 1 == getDimension());
   }

   std::set<int> MeshBase::getAttributes() const
   {
      return std::set<int>(
            getHandle().attributes.begin(), getHandle().attributes.end());
   }

   std::set<int> MeshBase::getBoundaryAttributes() const
   {
      return std::set<int>(
            getHandle().bdr_attributes.begin(), getHandle().bdr_attributes.end());
   }

   void Mesh<Context::Serial>::save(
         const boost::filesystem::path& filename, IO::FileFormat fmt, size_t precision) const
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

   MeshBase& MeshBase::update()
   {
      getHandle().SetAttributes();
      return *this;
   }

   std::deque<std::set<int>> MeshBase::ccl(
         std::function<bool(const Element&, const Element&)> p) const
   {
      std::set<int> visited;
      std::deque<int> searchQueue;
      std::deque<std::set<int>> res;

      // Perform the labelling
      assert(false);
      // for (int i = 0; i < count<Element>(); i++)
      // {
      //    if (!visited.count(i))
      //    {
      //       res.push_back({});
      //       searchQueue.push_back(i);
      //       while (searchQueue.size() > 0)
      //       {
      //          int el = searchQueue.back();
      //          searchQueue.pop_back();
      //          auto result = visited.insert(el);
      //          bool inserted = result.second;
      //          if (inserted)
      //          {
      //             res.back().insert(el);
      //             for (int n : get<Element>(el).adjacent())
      //             {
      //                if (p(get<Element>(el), get<Element>(n)))
      //                {
      //                   searchQueue.push_back(n);
      //                }
      //             }
      //          }
      //       }
      //    }
      // }
      return res;
   }

#ifdef RODIN_USE_MPI
   Mesh<Context::Parallel>
   Mesh<Context::Serial>::parallelize(boost::mpi::communicator comm)
   {
      return Mesh<Context::Parallel>(comm, *this);
   }
#endif

   // ---- Mesh<Serial> ------------------------------------------------------
   Mesh<Context::Serial>::Mesh(mfem::Mesh&& mesh)
      : m_mesh(std::move(mesh))
   {}

   Mesh<Context::Serial>::Mesh(const Mesh& other)
      : m_mesh(other.m_mesh)
   {}

   size_t Mesh<Context::Serial>::getCount(size_t dimension) const
   {
      if (dimension == getDimension())
      {
         return m_mesh.GetNE();
      }
      else if (dimension == getDimension() - 1)
      {
         return m_mesh.GetNumFaces();
      }
      else if (dimension == 0)
      {
         return m_mesh.GetNV();
      }
      else
      {
         assert(false);
      }
      return 0;
   }

   size_t Mesh<Context::Serial>::getElementCount() const
   {
      return getCount(getDimension());
   }

   BoundaryIterator Mesh<Context::Serial>::getBoundary() const
   {
      size_t dimension = getDimension() - 1;
      std::vector<Index> indices;
      indices.reserve(getHandle().GetNBE());
      for (int i = 0; i < getHandle().GetNBE(); i++)
      {
         int idx = getHandle().GetBdrFace(i);
         if (!getHandle().FaceIsInterior(idx))
            indices.push_back(idx);
      }
      return BoundaryIterator({dimension, *this, std::move(indices)});
   }

   InterfaceIterator Mesh<Context::Serial>::getInterface() const
   {
      size_t dimension = getDimension() - 1;
      std::vector<Index> indices;
      indices.reserve(getHandle().GetNumFaces());
      for (int idx = 0; idx < getHandle().GetNumFaces(); idx++)
      {
         if (getHandle().FaceIsInterior(idx))
            indices.push_back(idx);
      }
      return InterfaceIterator({dimension, *this, std::move(indices)});
   }

   ElementIterator Mesh<Context::Serial>::getElement(Index idx) const
   {
      size_t dimension = getDimension();
      std::vector<Index> indices(getCount(dimension));
      std::iota(indices.begin(), indices.end(), idx);
      return ElementIterator({dimension, *this, std::move(indices)});
   }

   FaceIterator Mesh<Context::Serial>::getFace(Index idx) const
   {
      size_t dimension = getDimension() - 1;
      std::vector<Index> indices(getCount(dimension));
      std::iota(indices.begin(), indices.end(), idx);
      return FaceIterator({dimension, *this, std::move(indices)});
   }

   VertexIterator Mesh<Context::Serial>::getVertex(Index idx) const
   {
      assert(false);
   }

   SimplexIterator Mesh<Context::Serial>::getSimplex(size_t dimension, Index idx) const
   {
      std::vector<Index> indices(getCount(dimension));
      std::iota(indices.begin(), indices.end(), idx);
      return SimplexIterator({dimension, *this, std::move(indices)});
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
         if (m_f2b.count(index))
         {
            return getHandle().GetBdrAttribute(m_f2b.at(index));
         }
         else
         {
            return RODIN_DEFAULT_SIMPLEX_ATTRIBUTE;
         }
      }
      else if (dimension == 0)
      {
         assert(false);
         return RODIN_DEFAULT_SIMPLEX_ATTRIBUTE;
      }
      else
      {
         assert(false);
         return RODIN_DEFAULT_SIMPLEX_ATTRIBUTE;
      }
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

   SubMesh<Context::Serial> Mesh<Context::Serial>::keep(int attr)
   {
      return keep(std::set<int>{attr});
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::keep(const std::set<int>& attrs)
   {
      // assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!
      assert(false);
      SubMesh<Context::Serial> res(*this);
      res.initialize(getDimension(), getSpaceDimension());
      std::vector<Index> indices;
      for (auto it = getElement(); !it.end(); ++it)
      {
         if (attrs.count(it->getAttribute()))
            indices.push_back(it->getIndex());
      }
      for (auto bit = getBoundary(); !bit.end(); ++bit)
      {
         for (auto eit = bit->getIncident(); !eit.end(); ++eit)
         {
            if (attrs.count(eit->getAttribute()))
               indices.push_back(eit->getIndex());
         }
      }
      res.include(getDimension(), indices);
      res.finalize();
      return res;
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::skin() const
   {
      assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!
      assert(false);
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::trim(int attr)
   {
      return trim(std::set<int>{attr});
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::trim(const std::set<int>& attrs)
   {
      std::set<int> complement = getAttributes();
      for (const auto& a : attrs)
         complement.erase(a);
      return keep(complement);
   }

   mfem::Mesh& Mesh<Context::Serial>::getHandle()
   {
      return m_mesh;
   }

   const mfem::Mesh& Mesh<Context::Serial>::getHandle() const
   {
      return m_mesh;
   }

   Mesh<Context::Serial>&
   Mesh<Context::Serial>::initialize(size_t dim, size_t sdim)
   {
      m_mesh = mfem::Mesh(dim, 0, 0, 0, sdim);
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::vertex(const std::vector<double>& x)
   {
      if (x.size() != getSpaceDimension())
      {
         Alert::Exception()
            << "Vertex dimension is different from space dimension"
            << " (" << x.size() << " != " << getSpaceDimension() << ")"
            << Alert::Raise;
      }
      getHandle().AddVertex(x.data());
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::element(
         Type geom,
         const std::vector<int>& vs, Attribute attr)
   {
      mfem::Element* el = getHandle().NewElement(static_cast<int>(geom));
      el->SetVertices(vs.data());
      el->SetAttribute(attr);
      getHandle().AddElement(el);
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::face(
         Type geom,
         const std::vector<int>& vs, Attribute attr)
   {
      mfem::Element* el = getHandle().NewElement(static_cast<int>(geom));
      el->SetVertices(vs.data());
      el->SetAttribute(attr);
      getHandle().AddBdrElement(el);
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::finalize()
   {
      getHandle().FinalizeTopology();
      getHandle().Finalize(false, true);
      for (int i = 0; i < getHandle().GetNBE(); i++)
         m_f2b[getHandle().GetBdrElementEdgeIndex(i)] = i;
      return *this;
   }
}

