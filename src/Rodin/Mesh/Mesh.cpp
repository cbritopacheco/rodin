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

namespace Rodin
{
   // ---- MeshBase ----------------------------------------------------------
   int MeshBase::getSpaceDimension() const
   {
      return getHandle().SpaceDimension();
   }

   int MeshBase::getDimension() const
   {
      return getHandle().Dimension();
   }

   bool MeshBase::isSurface() const
   {
      return (getSpaceDimension() - 1 == getDimension());
   }

   void MeshBase::refine()
   {
      getHandle().UniformRefinement();
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
      assert(u.getFiniteElementSpace().getVectorDimension() == getDimension());
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
      for (int i = 0; i < count<Element>(); i++)
         totalVolume += getHandle().GetElementVolume(i);
      return totalVolume;
   }

   double MeshBase::getVolume(int attr)
   {
      double totalVolume = 0;
      for (int i = 0; i < count<Element>(); i++)
         totalVolume += getHandle().GetElementVolume(i) * (getHandle().GetAttribute(i) == attr);
      return totalVolume;
   }

   double MeshBase::getBoundaryElementArea(int i)
   {
      mfem::ElementTransformation *et = getHandle().GetBdrElementTransformation(i);
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(
            getHandle().GetBdrElementBaseGeometry(i), et->OrderJ());
      double area = 0.0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(j);
         et->SetIntPoint(&ip);
         area += ip.weight * et->Weight();
      }
      return area;
   }

   double MeshBase::getPerimeter()
   {
      double totalArea = 0;
      for (int i = 0; i < count<Element>(); i++)
         totalArea += getBoundaryElementArea(i);
      return totalArea;
   }

   double MeshBase::getPerimeter(int attr)
   {
      double totalVolume = 0;
      for (int i = 0; i < count<BoundaryElement>(); i++)
         totalVolume += getBoundaryElementArea(i) * (getHandle().GetBdrAttribute(i) == attr);
      return totalVolume;
   }

   std::set<int> MeshBase::where(std::function<bool(const Element&)> p) const
   {
      std::set<int> res;
      for (int i = 0; i < count<Element>(); i++)
         if (p(get<Element>(i)))
            res.insert(i);
      return res;
   }

   MeshBase& MeshBase::edit(std::function<void(ElementView)> f)
   {
      for (int i = 0; i < count<Element>(); i++)
         f(get<Element>(i));
      return *this;
   }

   MeshBase& MeshBase::edit(std::function<void(BoundaryElementView)> f)
   {
      for (int i = 0; i < count<BoundaryElement>(); i++)
         f(get<BoundaryElement>(i));
      return *this;
   }

   MeshBase& MeshBase::edit(std::function<void(ElementView)> f, const std::set<int>& elements)
   {
      for (auto el : elements)
      {
         assert(el >= 0);
         assert(el < count<Element>());

         f(get<Element>(el));
      }
      return *this;
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
      for (int i = 0; i < count<Element>(); i++)
      {
         if (!visited.count(i))
         {
            res.push_back({});
            searchQueue.push_back(i);
            while (searchQueue.size() > 0)
            {
               int el = searchQueue.back();
               searchQueue.pop_back();
               auto result = visited.insert(el);
               bool inserted = result.second;
               if (inserted)
               {
                  res.back().insert(el);
                  for (int n : get<Element>(el).adjacent())
                  {
                     if (p(get<Element>(el), get<Element>(n)))
                     {
                        searchQueue.push_back(n);
                     }
                  }
               }
            }
         }
      }
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
         }
      }

      return *this;
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::extract(const std::set<int>& elements)
   {
      SubMesh<Context::Serial> res(*this);
      res.initialize(getDimension(), getSpaceDimension(), elements.size());
      for (int el : elements)
         res.add(get<Element>(el));
      res.finalize();
      return res;
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::keep(int attr)
   {
      return keep(std::set<int>{attr});
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::keep(const std::set<int>& attrs)
   {
      assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!

      SubMesh<Context::Serial> res(*this);
      res.initialize(getDimension(), getSpaceDimension());

      // Add elements with matching attribute
      for (int i = 0; i < count<Element>(); i++)
      {
         const auto& el = get<Element>(i);
         if (attrs.count(el.getAttribute()))
            res.add(el);
      }

      // Add the boundary elements
      for (int i = 0; i < count<BoundaryElement>(); i++)
      {
         const auto& be = get<BoundaryElement>(i);
         const auto& elems = be.elements();
         for (const auto& el : elems)
         {
            if (attrs.count(get<Element>(el).getAttribute()))
            {
               res.add(be);
               break;
            }
         }
      }
      res.finalize();
      return res;
   }

   SubMesh<Context::Serial> Mesh<Context::Serial>::skin()
   {
      assert(!getHandle().GetNodes()); // Curved mesh or discontinuous mesh not handled yet!

      SubMesh<Context::Serial> res(*this);
      res.initialize(getSpaceDimension() - 1, getSpaceDimension());
      for (int i = 0; i < count<BoundaryElement>(); i++)
         res.add(get<BoundaryElement>(i));
      res.finalize();
      return res;
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
      for (int i = 0; i < count<Face>(); i++)
      {
         const auto& fc = get<Face>(i);
         const auto& elems = fc.elements();
         if (elems.size() == 2)
         {
            std::set<int> k{
               get<Element>(*elems.begin()).getAttribute(),
               get<Element>(*std::next(elems.begin())).getAttribute()
            };
            auto it = boundaries.find(k);
            if (it != boundaries.end())
            {
               mfem::Element* be =
                  getHandle().NewElement(fc.getHandle().GetGeometryType());
               be->SetVertices(fc.getHandle().GetVertices());
               be->SetAttribute(it->second);
               getHandle().AddBdrElement(be);
            }
         }
      }
      getHandle().FinalizeTopology();
      return *this;
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
   Mesh<Context::Serial>::initialize(int dim, int sdim, int nv)
   {
      m_mesh = mfem::Mesh(dim, nv, 0, 0, sdim);
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
         Geometry geom,
         const std::vector<int>& vs, std::optional<int> attr)
   {
      mfem::Element* el = getHandle().NewElement(static_cast<int>(geom));
      el->SetVertices(vs.data());
      if (attr)
         el->SetAttribute(*attr);
      getHandle().AddElement(el);
      return *this;
   }

   Mesh<Context::Serial>& Mesh<Context::Serial>::finalize()
   {
      getHandle().Finalize(false, true);
      return *this;
   }

   // ---- Mesh<Parallel> ----------------------------------------------------
}

