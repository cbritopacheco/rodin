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
   ElementView MeshBase::getElement(int i)
   {
      return ElementView(*this, getHandle().GetElement(i), i);
   }

   Element MeshBase::getElement(int i) const
   {
      return Element(*this, getHandle().GetElement(i), i);
   }

   BoundaryElementView MeshBase::getBoundaryElement(int i)
   {
      return BoundaryElementView(*this, getHandle().GetBdrElement(i), i);
   }

   BoundaryElement MeshBase::getBoundaryElement(int i) const
   {
      return BoundaryElement(*this, getHandle().GetBdrElement(i), i);
   }

   int MeshBase::getElementCount() const
   {
      return getHandle().GetNE();
   }

   int MeshBase::getBoundaryElementCount() const
   {
      return getHandle().GetNBE();
   }

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

   void Mesh<Traits::Serial>::save(
         const boost::filesystem::path& filename, IO::MeshFormat fmt, int precision) const
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
         case IO::MeshFormat::MFEM:
         {
            IO::MeshPrinter<IO::MeshFormat::MFEM, Traits::Serial> printer(*this);
            printer.print(ofs);
            break;
         }
         case IO::MeshFormat::GMSH:
         {
            IO::MeshPrinter<IO::MeshFormat::GMSH, Traits::Serial> printer(*this);
            printer.print(ofs);
            break;
         }
         case IO::MeshFormat::MEDIT:
         {
            IO::MeshPrinter<IO::MeshFormat::MEDIT, Traits::Serial> printer(*this);
            printer.print(ofs);
            break;
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
      for (int i = 0; i < getElementCount(); i++)
         totalVolume += getHandle().GetElementVolume(i);
      return totalVolume;
   }

   double MeshBase::getVolume(int attr)
   {
      double totalVolume = 0;
      for (int i = 0; i < getElementCount(); i++)
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
      for (int i = 0; i < getElementCount(); i++)
         totalArea += getBoundaryElementArea(i);
      return totalArea;
   }

   double MeshBase::getPerimeter(int attr)
   {
      double totalVolume = 0;
      for (int i = 0; i < getBoundaryElementCount(); i++)
         totalVolume += getBoundaryElementArea(i) * (getHandle().GetBdrAttribute(i) == attr);
      return totalVolume;
   }

   std::set<int> MeshBase::where(std::function<bool(const Element&)> p) const
   {
      std::set<int> res;
      for (int i = 0; i < getElementCount(); i++)
         if (p(getElement(i)))
            res.insert(i);
      return res;
   }

   MeshBase& MeshBase::edit(std::function<void(ElementView)> f)
   {
      for (int i = 0; i < getElementCount(); i++)
         f(getElement(i));
      return *this;
   }

   MeshBase& MeshBase::edit(std::function<void(ElementView)> f, const std::set<int>& elements)
   {
      for (auto el : elements)
      {
         assert(el >= 0);
         assert(el < getElementCount());

         f(getElement(el));
      }
      return *this;
   }

   std::deque<std::set<int>> MeshBase::ccl(
         std::function<bool(const Element&, const Element&)> p) const
   {
      std::set<int> visited;
      std::deque<int> searchQueue;
      std::deque<std::set<int>> res;

      // Perform the labelling
      for (int i = 0; i < getElementCount(); i++)
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
                  for (int n : getElement(el).adjacent())
                  {
                     if (p(getElement(el), getElement(n)))
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
   Mesh<Traits::Parallel>
   Mesh<Traits::Serial>::parallelize(boost::mpi::communicator comm)
   {
      return Mesh<Traits::Parallel>(comm, *this);
   }
#endif

   // ---- Mesh<Serial> ------------------------------------------------------
   Mesh<Traits::Serial>::Mesh(mfem::Mesh&& mesh)
      : m_mesh(std::move(mesh))
   {}

   Mesh<Traits::Serial>::Mesh(const Mesh& other)
      : m_mesh(other.m_mesh)
   {}

   Mesh<Traits::Serial>&
   Mesh<Traits::Serial>::load(const boost::filesystem::path& filename, IO::MeshFormat fmt)
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
         case IO::MeshFormat::MFEM:
         {
            IO::MeshLoader<IO::MeshFormat::MFEM, Traits::Serial> loader(*this);
            loader.load(input);
            break;
         }
         case IO::MeshFormat::GMSH:
         {
            IO::MeshLoader<IO::MeshFormat::GMSH, Traits::Serial> loader(*this);
            loader.load(input);
            break;
         }
         case IO::MeshFormat::MEDIT:
         {
            IO::MeshLoader<IO::MeshFormat::MEDIT, Traits::Serial> loader(*this);
            loader.load(input);
            break;
         }
      }

      return *this;
   }

   SubMesh<Traits::Serial> Mesh<Traits::Serial>::trim(int attr, int bdrLabel)
   {
      return trim(std::set<int>{attr}, bdrLabel);
   }

   SubMesh<Traits::Serial> Mesh<Traits::Serial>::trim(const std::set<int>& attrs, int bdrLabel)
   {
      assert(attrs.size() > 0);

      // Count the number of elements in the trimmed mesh
      int ne = 0;
      for (int i = 0; i < getElementCount(); i++)
      {
         int elemAttr = getHandle().GetElement(i)->GetAttribute();
         ne += !attrs.count(elemAttr); // Count is always 0 or 1
      }

      // Count the number of boundary elements in the trimmed mesh
      int nbe = 0;
      for (int i = 0; i < getHandle().GetNumFaces(); i++)
      {
         int e1 = -1, e2 = -1;
         getHandle().GetFaceElements(i, &e1, &e2);

         int a1 = 0, a2 = 0;
         if (e1 >= 0)
            a1 = getHandle().GetElement(e1)->GetAttribute();
         if (e2 >= 0)
            a2 = getHandle().GetElement(e2)->GetAttribute();

         if (a1 == 0 || a2 == 0)
         {
            nbe += (a1 == 0) && !attrs.count(a2);
            nbe += (a2 == 0 && !attrs.count(a1));
         }
         else
         {
            nbe += attrs.count(a1) && !attrs.count(a2);
            nbe += !attrs.count(a1) && attrs.count(a2);
         }
      }

      mfem::Mesh trimmed(getHandle().Dimension(), getHandle().GetNV(),
            ne, nbe, getHandle().SpaceDimension());

      // Copy vertices
      for (int i = 0; i < getHandle().GetNV(); i++)
         trimmed.AddVertex(getHandle().GetVertex(i));

      // Copy elements
      for (int i = 0; i < getElementCount(); i++)
      {
         mfem::Element* el = getHandle().GetElement(i);
         int elemAttr = el->GetAttribute();
         if (!attrs.count(elemAttr))
         {
            mfem::Element* nel = getHandle().NewElement(el->GetGeometryType());
            nel->SetAttribute(elemAttr);
            nel->SetVertices(el->GetVertices());
            trimmed.AddElement(nel);
         }
      }

      // Copy selected boundary elements
      for (int i = 0; i < getBoundaryElementCount(); i++)
      {
         int e, info;
         getHandle().GetBdrElementAdjacentElement(i, e, info);

         int elemAttr = getHandle().GetElement(e)->GetAttribute();
         if (!attrs.count(elemAttr))
         {
            mfem::Element* nbel = getHandle().GetBdrElement(i)->Duplicate(&trimmed);
            trimmed.AddBdrElement(nbel);
         }
      }

      // Create new boundary elements
      for (int i = 0; i < getHandle().GetNumFaces(); i++)
      {
         int e1 = -1, e2 = -1;
         getHandle().GetFaceElements(i, &e1, &e2);

         int i1 = -1, i2 = -1;
         getHandle().GetFaceInfos(i, &i1, &i2);

         int a1 = 0, a2 = 0;
         if (e1 >= 0)
            a1 = getHandle().GetElement(e1)->GetAttribute();
         if (e2 >= 0)
            a2 = getHandle().GetElement(e2)->GetAttribute();

         if (a1 != 0 && a2 != 0)
         {
            if (attrs.count(a1) && !attrs.count(a2))
            {
               mfem::Element* bel;
               if (getHandle().Dimension() == 1)
                  bel = new mfem::Point(&i);
               else
                  bel = getHandle().GetFace(i)->Duplicate(&trimmed);
               bel->SetAttribute(bdrLabel);
               trimmed.AddBdrElement(bel);
            }
            else if (!attrs.count(a1) && attrs.count(a2))
            {
               mfem::Element* bel;
               if (getHandle().Dimension() == 1)
                  bel = new mfem::Point(&i);
               else
                  bel = getHandle().GetFace(i)->Duplicate(&trimmed);
               bel->SetAttribute(bdrLabel);
               trimmed.AddBdrElement(bel);
            }
         }
      }

      trimmed.FinalizeTopology();
      trimmed.Finalize();

      /* Get map of vertices from submesh to mesh
       *
       * The following code is taken from RemoveUnusedVertices() method and
       * slightly modified to be clearer. It computes the submesh to mesh map
       * of the vertex numberings.
       */
      mfem::Array<int> vertexUsage(trimmed.GetNV());
      vertexUsage = 0;
      for (int i = 0; i < trimmed.GetNE(); i++)
      {
         mfem::Element* el = trimmed.GetElement(i);
         const int* v = el->GetVertices();
         for (int j = 0; j < el->GetNVertices(); j++)
         {
            assert(v[j] < trimmed.GetNV());
            vertexUsage[v[j]] = 1;
         }
      }
      for (int i = 0; i < trimmed.GetNBE(); i++)
      {
         mfem::Element* el = trimmed.GetBdrElement(i);
         const int* v = el->GetVertices();
         for (int j = 0; j < el->GetNVertices(); j++)
         {
            assert(v[j] < trimmed.GetNV());
            vertexUsage[v[j]] = 1;
         }
      }
      int vertexIdx = 0;
      boost::bimap<int, int> s2pv; // Submesh to mesh node map
      for (int i = 0; i < vertexUsage.Size(); i++)
      {
         if (vertexUsage[i])
            s2pv.insert({vertexIdx++, i});
      }

      trimmed.RemoveUnusedVertices();

      // Check for curved or discontinuous mesh
      if (getHandle().GetNodes())
      {
         // Extract Nodes GridFunction and determine its type
         const mfem::GridFunction* Nodes = getHandle().GetNodes();
         const mfem::FiniteElementSpace* fes = Nodes->FESpace();

         mfem::Ordering::Type ordering = fes->GetOrdering();
         int order = fes->FEColl()->GetOrder();
         int sdim = getHandle().SpaceDimension();
         bool discont =
            dynamic_cast<const mfem::L2_FECollection*>(fes->FEColl()) != NULL;

         // Set curvature of the same type as original mesh
         trimmed.SetCurvature(order, discont, sdim, ordering);

         const mfem::FiniteElementSpace* trimmedFES = trimmed.GetNodalFESpace();
         mfem::GridFunction* trimmedNodes = trimmed.GetNodes();

         mfem::Array<int> vdofs;
         mfem::Array<int> trimmedVdofs;
         mfem::Vector locVec;

         // Copy nodes to trimmed mesh
         int te = 0;
         for (int e = 0; e < getElementCount(); e++)
         {
            mfem::Element* el = getHandle().GetElement(e);
            int elemAttr = el->GetAttribute();
            if (!attrs.count(elemAttr))
            {
               fes->GetElementVDofs(e, vdofs);
               Nodes->GetSubVector(vdofs, locVec);

               trimmedFES->GetElementVDofs(te, trimmedVdofs);
               trimmedNodes->SetSubVector(trimmedVdofs, locVec);
               te++;
            }
         }
      }
      return SubMesh<Traits::Serial>(
            std::move(trimmed)).setParent(*this).setVertexMap(std::move(s2pv));
   }

   SubMesh<Traits::Serial> Mesh<Traits::Serial>::skin(const std::map<int, int>& bdrAttr)
   {
      auto res = skin();
      if (bdrAttr.size() == 0)
         return res;

      for (int i = 0; i < res.getHandle().GetNumFaces(); i++)
      {
         int el1, el2;
         res.getHandle().GetFaceElements(i, &el1, &el2);

         int attr1 = res.getHandle().GetAttribute(el1);
         int attr2 = res.getHandle().GetAttribute(el2);

         if (attr1 != attr2)
         {
            int attr;
            if (bdrAttr.count(attr1))
               attr = bdrAttr.at(attr1);
            else if (bdrAttr.count(attr2))
               attr = bdrAttr.at(attr2);
            else
               continue;

            mfem::Element* fc = res.getHandle().GetFace(i)->Duplicate(&res.getHandle());
            fc->SetAttribute(attr);
            res.getHandle().AddBdrElement(fc);
         }
      }
      res.getHandle().SetAttributes();

      return res;
   }

   SubMesh<Traits::Serial> Mesh<Traits::Serial>::skin()
   {
      assert(getDimension() >= 2);

      // Determine mapping from vertex to boundary vertex
      std::set<int> bdrVertices;
      for (int i = 0; i < getBoundaryElementCount(); i++)
      {
         mfem::Element* el = getHandle().GetBdrElement(i);
         int *v = el->GetVertices();
         int nv = el->GetNVertices();
         for (int j = 0; j < nv; j++)
            bdrVertices.insert(v[j]);
      }

      mfem::Mesh res(getDimension() - 1, bdrVertices.size(),
            getBoundaryElementCount(), 0, getSpaceDimension());

      // Copy vertices to the boundary mesh
      int vertexIdx = 0;
      boost::bimap<int, int> s2pv; // Submesh to mesh node map
      for (const auto& v : bdrVertices)
      {
         s2pv.insert({vertexIdx, v});
         vertexIdx++;
         double *c = getHandle().GetVertex(v);
         res.AddVertex(c);
      }

      // Copy elements to the boundary mesh
      for (int i = 0; i < getBoundaryElementCount(); i++)
      {
         mfem::Element *el = getHandle().GetBdrElement(i);
         int *v = el->GetVertices();
         int nv = el->GetNVertices();

         std::vector<int> bv(nv);
         for (int j = 0; j < nv; j++)
            bv[j] = s2pv.right.at(v[j]);

         switch (el->GetGeometryType())
         {
            case mfem::Geometry::SEGMENT:
               res.AddSegment(bv.data(), el->GetAttribute());
               break;
            case mfem::Geometry::TRIANGLE:
               res.AddTriangle(bv.data(), el->GetAttribute());
               break;
            case mfem::Geometry::SQUARE:
               res.AddQuad(bv.data(), el->GetAttribute());
               break;
            default:
               break; // This should not happen
         }
      }
      res.FinalizeTopology();

      // Copy GridFunction describing nodes if present
      if (getHandle().GetNodes())
      {
         mfem::FiniteElementSpace* fes = getHandle().GetNodes()->FESpace();
         const mfem::FiniteElementCollection* fec = fes->FEColl();
         if (dynamic_cast<const mfem::H1_FECollection*>(fec))
         {
            mfem::FiniteElementCollection *fec_copy =
               mfem::FiniteElementCollection::New(fec->Name());
            mfem::FiniteElementSpace* fes_copy =
               new mfem::FiniteElementSpace(*fes, &res, fec_copy);
            mfem::GridFunction *bdr_nodes = new mfem::GridFunction(fes_copy);
            bdr_nodes->MakeOwner(fec_copy);

            res.NewNodes(*bdr_nodes, true);

            mfem::Array<int> vdofs;
            mfem::Array<int> bvdofs;
            mfem::Vector v;
            for (int i = 0; i < getBoundaryElementCount(); i++)
            {
               fes->GetBdrElementVDofs(i, vdofs);
               getHandle().GetNodes()->GetSubVector(vdofs, v);

               fes_copy->GetElementVDofs(i, bvdofs);
               bdr_nodes->SetSubVector(bvdofs, v);
            }
         }
         else
         {
            Alert::Exception("Discontinuous node space not yet supported.").raise();
         }
      }
      return SubMesh<Traits::Serial>(
            std::move(res)).setVertexMap(std::move(s2pv)).setParent(*this);
   }

   mfem::Mesh& Mesh<Traits::Serial>::getHandle()
   {
      return m_mesh;
   }

   const mfem::Mesh& Mesh<Traits::Serial>::getHandle() const
   {
      return m_mesh;
   }

   // ---- Mesh<Parallel> ----------------------------------------------------
}

