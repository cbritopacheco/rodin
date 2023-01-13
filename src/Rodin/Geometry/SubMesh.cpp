#include "SubMesh.h"

#include "Element.h"

namespace Rodin::Geometry
{
   SubMesh<Context::Serial>::SubMesh(const MeshBase& parent)
      : m_parent(parent)
   {}

   SubMesh<Context::Serial>::SubMesh(const SubMesh& other)
      :  Mesh(other),
         m_parent(other.m_parent),
         m_s2pv(other.m_s2pv)
   {}

   const MeshBase& SubMesh<Context::Serial>::getParent() const
   {
      return m_parent;
   }

   SubMesh<Context::Serial>&
   SubMesh<Context::Serial>::include(size_t dim, const std::set<Index>& simplices)
   {
      if (dim == getDimension())
      {
         for (const auto& simplex : simplices)
         {
            auto el = getParent().getElement(simplex);

            // Add element vertices to the resulting mesh
            mfem::Array<int> pv;
            el->getHandle().GetVertices(pv);

            mfem::Array<int> sv(pv.Size());
            for (int i = 0; i < sv.Size(); i++)
            {
               int pvid = pv[i];
               if (m_s2pv.right.count(pvid) == 0) // Only add vertex if it is not in the map
               {
                  sv[i] = getHandle().AddVertex(getParent().getHandle().GetVertex(pvid));
                  m_s2pv.insert({static_cast<size_t>(sv[i]), static_cast<size_t>(pvid)});
               }
               else // Else get the id of the vertex in the submesh
               {
                  sv[i] = m_s2pv.right.at(pvid);
               }
            }

            // Add element with the new vertex ordering
            mfem::Element* newEl =
               getHandle().NewElement(el->getHandle().GetGeometryType());
            newEl->SetVertices(sv);
            newEl->SetAttribute(el->getAttribute());
            size_t seid = getHandle().AddElement(newEl);
            m_s2pe.insert({seid, el->getIndex()});
         }

         for (int i = 0; i < getParent().getHandle().GetNBE(); i++)
         {
            Index faceIdx = getParent().getHandle().GetBdrFace(i);
            Attribute attr = getParent().getFaceAttribute(faceIdx);
            int el1 = -1, el2 = -1;
            getParent().getHandle().GetFaceElements(faceIdx, &el1, &el2);
            if (simplices.count(el1) || simplices.count(el2))
            {
               assert(el1 >= 0 || el2 >= 0);
               mfem::Array<int> vs;
               getParent().getHandle().GetFaceVertices(faceIdx, vs);
               for (int j = 0; j < vs.Size(); j++)
                  vs[j] = m_s2pv.right.at(vs[j]);
               mfem::Element* newEl =
                  getHandle().NewElement(getParent().getHandle().GetFaceGeometry(faceIdx));
               newEl->SetVertices(vs);
               newEl->SetAttribute(attr);
               getHandle().AddBdrElement(newEl);
            }
         }
      }
      else if (dim == getDimension() - 1)
      {
         assert(false);
      }
      else
      {
         assert(false);
      }

      return *this;
   }

   const boost::bimap<size_t, size_t>& SubMesh<Context::Serial>::getVertexMap() const
   {
      return m_s2pv;
   }

   const boost::bimap<size_t, size_t>& SubMesh<Context::Serial>::getElementMap() const
   {
      return m_s2pe;
   }
}
