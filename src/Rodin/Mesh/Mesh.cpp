/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Mesh.h"
#include "SubMesh.h"

namespace Rodin
{
   Mesh Mesh::load(const std::filesystem::path& filename)
   {
      return Mesh(mfem::Mesh::LoadFromFile(filename.c_str()));
   }

   void Mesh::save(const std::filesystem::path& filename)
   {
      getHandle().Save(filename.c_str());
   }

   Mesh::Mesh(mfem::Mesh&& mesh)
      : m_mesh(std::move(mesh))
   {}

   Mesh::Mesh(const Mesh& other)
      : m_mesh(other.m_mesh)
   {}

   Mesh& Mesh::displace(const Variational::GridFunctionBase& u)
   {
      assert(u.getFiniteElementSpace().getRangeDimension() == getDimension());
      m_mesh.MoveNodes(u.getHandle());
      return *this;
   }

   double
   Mesh::getMaximumDisplacement(const Variational::GridFunctionBase& u)
   {
      double res;
      m_mesh.CheckDisplacements(u.getHandle(), res);
      return res;
   }

   int Mesh::getDimension() const
   {
      return m_mesh.Dimension();
   }

   double Mesh::getVolume()
   {
      double totalVolume = 0;
      for (int i = 0; i < m_mesh.GetNE(); i++)
         totalVolume += m_mesh.GetElementVolume(i);
      return totalVolume;
   }

   double Mesh::getVolume(int attr)
   {
      double totalVolume = 0;
      for (int i = 0; i < m_mesh.GetNE(); i++)
      {
         if (m_mesh.GetAttribute(i) == attr)
            totalVolume += m_mesh.GetElementVolume(i);
      }
      return totalVolume;
   }

   SubMesh Mesh::trim(int attr, int bdrLabel)
   {
      return trim(std::set<int>{attr}, bdrLabel);
   }

   SubMesh Mesh::trim(const std::set<int>& attrs, int bdrLabel)
   {
      assert(attrs.size() > 0);

      // Count the number of elements in the trimmed mesh
      int ne = 0;
      for (int i = 0; i < getHandle().GetNE(); i++)
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
      for (int i = 0; i < getHandle().GetNE(); i++)
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
      for (int i = 0; i < getHandle().GetNBE(); i++)
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
      std::map<int, int> s2pv; // Submesh to mesh node map
      for (int i = 0; i < vertexUsage.Size(); i++)
      {
         if (vertexUsage[i])
            s2pv[vertexIdx++] = i;
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
         for (int e = 0; e < getHandle().GetNE(); e++)
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
      return SubMesh(std::move(trimmed)).setParent(*this).setVertexMap(std::move(s2pv));
   }

   mfem::Mesh& Mesh::getHandle()
   {
      return m_mesh;
   }

   const mfem::Mesh& Mesh::getHandle() const
   {
      return m_mesh;
   }
}

