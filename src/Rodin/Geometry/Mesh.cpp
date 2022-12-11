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
         const boost::filesystem::path& filename, IO::FileFormat fmt, int precision) const
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
      for (auto it = begin(getDimension()); !it.end(); it++)
         totalVolume += it->getVolume();
      return totalVolume;
   }

   double MeshBase::getVolume(Attribute attr)
   {
      double totalVolume = 0;
      for (auto it = getElement(getDimension()); !it.end(); it++)
      {
         if (it->getAttribute() == attr)
            totalVolume += it->getVolume();
      }
      return totalVolume;
   }

   double MeshBase::getPerimeter()
   {
      double totalVolume = 0;
      for (auto it = getBoundary(); !it.end(); it++)
         totalVolume += it->getVolume();
      return totalVolume;
   }

   double MeshBase::getPerimeter(Attribute attr)
   {
      double totalVolume = 0;
      for (auto it = getBoundary(); !it.end(); it++)
      {
         if (it->getAttribute() == attr)
            totalVolume += it->getVolume();
      }
      return totalVolume;
   }

   // std::set<int> MeshBase::where(std::function<bool(const Element&)> p) const
   // {
   //    std::set<int> res;
   //    for (int i = 0; i < count<Element>(); i++)
   //       if (p(get<Element>(i)))
   //          res.insert(i);
   //    return res;
   // }

   // MeshBase& MeshBase::edit(std::function<void(ElementView)> f)
   // {
   //    for (int i = 0; i < count<Element>(); i++)
   //       f(get<Element>(i));
   //    return *this;
   // }

   // MeshBase& MeshBase::edit(std::function<void(BoundaryView)> f)
   // {
   //    for (int i = 0; i < count<Boundary>(); i++)
   //       f(get<Boundary>(i));
   //    return *this;
   // }

   // MeshBase& MeshBase::edit(std::function<void(ElementView)> f, const std::set<int>& elements)
   // {
   //    for (auto el : elements)
   //    {
   //       assert(el >= 0);
   //       assert(el < count<Element>());

   //       f(get<Element>(el));
   //    }
   //    return *this;
   // }

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

   MeshSimplexIterator
   Mesh<Context::Serial>::begin(size_t dimension) noexcept
   {
      assert(false);
      // return SimplexIterator(*this, dimension, 0);
   }

   const MeshSimplexIterator
   Mesh<Context::Serial>::begin(size_t dimension) const noexcept
   {
      assert(false);
      // return SimplexIterator(*this, dimension, 0);
   }

   const MeshSimplexIterator
   Mesh<Context::Serial>::cbegin(size_t dimension) const noexcept
   {
      assert(false);
      // return SimplexIterator(*this, dimension, 0);
   }

   MeshSimplexIterator
   Mesh<Context::Serial>::end(size_t dimension) noexcept
   {
      assert(false);
      // return SimplexIterator(*this, dimension, getCount(dimension), false);
   }

   const MeshSimplexIterator
   Mesh<Context::Serial>::end(size_t dimension) const noexcept
   {
      assert(false);
      // return SimplexIterator(*this, dimension, getCount(dimension), false);
   }

   const MeshSimplexIterator
   Mesh<Context::Serial>::cend(size_t dimension) const noexcept
   {
      assert(false);
      // return SimplexIterator(*this, dimension, getCount(dimension), false);
   }

   size_t Mesh<Context::Serial>::getCount(Region region) const
   {
      switch (region)
      {
         case Region::Domain:
         {
            return getHandle().GetNE();
         }
         case Region::Boundary:
         {
            assert(false);
            return 0;
         }
         case Region::Interface:
         {
            assert(false);
            return 0;
         }
      }
   }

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

   MeshBoundaryIterator Mesh<Context::Serial>::getBoundary() const
   {
      assert(false);
   }

   MeshInterfaceIterator Mesh<Context::Serial>::getInterface() const
   {
      assert(false);
   }

   MeshElementIterator Mesh<Context::Serial>::getElement(size_t idx) const
   {
      assert(false);
   }

   MeshFaceIterator Mesh<Context::Serial>::getFace(size_t idx) const
   {
      assert(false);
   }

   MeshVertexIterator Mesh<Context::Serial>::getVertex(size_t idx) const
   {
      assert(false);
   }

   MeshSimplexIterator Mesh<Context::Serial>::getSimplex(size_t idx, size_t dimension) const
   {
      assert(false);
   }

   bool Mesh<Context::Serial>::isInterface(Index faceIdx) const
   {
      return getHandle().FaceIsInterior(faceIdx);
   }

   bool Mesh<Context::Serial>::isBoundary(Index faceIdx) const
   {
      return !getHandle().FaceIsInterior(faceIdx);
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

   SubMesh<Context::Serial> Mesh<Context::Serial>::extract(const std::set<int>& elements)
   {
      // SubMesh<Context::Serial> res(*this);
      // res.initialize(getDimension(), getSpaceDimension(), elements.size());
      // for (int el : elements)
      //    res.add(get<Element>(el));
      // res.finalize();
      // return res;
      assert(false);
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::keep(int attr)
   {
      return keep(std::set<int>{attr});
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::keep(const std::set<int>& attrs)
   {
      // assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!

      // SubMesh<Context::Serial> res(*this);
      // res.initialize(getDimension(), getSpaceDimension());

      // // Add elements with matching attribute
      // for (int i = 0; i < count<Element>(); i++)
      // {
      //    const auto& el = get<Element>(i);
      //    if (attrs.count(el.getAttribute()))
      //       res.add(el);
      // }

      // // Add the boundary elements
      // for (int i = 0; i < count<Boundary>(); i++)
      // {
      //    const auto& be = get<Boundary>(i);
      //    const auto& elems = be.elements();
      //    for (const auto& el : elems)
      //    {
      //       if (attrs.count(get<Element>(el).getAttribute()))
      //       {
      //          res.add(be);
      //          break;
      //       }
      //    }
      // }
      // res.finalize();
      // return res;
      assert(false);
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::skin()
   {
      assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!

      assert(false);
      // SubMesh<Context::Serial> res(*this);
      // res.initialize(getSpaceDimension() - 1, getSpaceDimension());
      // for (int i = 0; i < count<Boundary>(); i++)
      //    res.add(get<Boundary>(i));
      // res.finalize();
      // return res;
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

   Mesh<Context::Serial>& Mesh<Context::Serial>::trace(
         const std::map<std::set<int>, int>& boundaries)
   {
      // for (int i = 0; i < count<Face>(); i++)
      // {
      //    const auto& fc = get<Face>(i);
      //    const auto& elems = fc.elements();
      //    if (elems.size() == 2)
      //    {
      //       std::set<int> k{
      //          get<Element>(*elems.begin()).getAttribute(),
      //          get<Element>(*std::next(elems.begin())).getAttribute()
      //       };
      //       auto it = boundaries.find(k);
      //       if (it != boundaries.end())
      //       {
      //          mfem::Element* be =
      //             getHandle().NewElement(fc.getHandle().GetGeometryType());
      //          be->SetVertices(fc.getHandle().GetVertices());
      //          be->SetAttribute(it->second);
      //          getHandle().AddBdrElement(be);
      //       }
      //    }
      // }
      // getHandle().FinalizeTopology();
      // return *this;
      assert(false);
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
   Mesh<Context::Serial>::initialize(int dim, int sdim)
   {
      m_mesh = mfem::Mesh(dim, 0, 0, 0, sdim);
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::vertex(const std::vector<double>& x)
   {
      if (static_cast<int>(x.size()) != getSpaceDimension())
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
         const std::vector<int>& vs, int attr)
   {
      mfem::Element* el = getHandle().NewElement(static_cast<int>(geom));
      el->SetVertices(vs.data());
      el->SetAttribute(attr);
      getHandle().AddElement(el);
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::boundary(
         Type geom,
         const std::vector<int>& vs, int attr)
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
      return *this;
   }

   // ---- Mesh<Parallel> ----------------------------------------------------
}

