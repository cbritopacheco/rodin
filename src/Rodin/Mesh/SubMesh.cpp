#include "SubMesh.h"

#include "Element.h"

namespace Rodin
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

   SubMesh<Context::Serial>& SubMesh<Context::Serial>::add(const Element& el)
   {
      assert(&getParent() == &el.getMesh());
      assert(getDimension() == getParent().getDimension());
      assert(getSpaceDimension() == getParent().getSpaceDimension());

      // Add element vertices to the resulting mesh
      mfem::Array<int> pv;
      el.getHandle().GetVertices(pv);

      mfem::Array<int> sv(pv.Size());
      for (int i = 0; i < sv.Size(); i++)
      {
         int pvid = pv[i];
         if (m_s2pv.right.count(pvid) == 0) // Only add vertex if it is not in the map
         {
            sv[i] = getHandle().AddVertex(getParent().getHandle().GetVertex(pvid));
            m_s2pv.insert({sv[i], pvid});
         }
         else // Else get the id of the vertex in the submesh
         {
            sv[i] = m_s2pv.right.at(pvid);
         }
      }

      // Add element with the new vertex ordering
      mfem::Element* newEl =
         getHandle().NewElement(el.getHandle().GetGeometryType());
      newEl->SetVertices(sv);
      newEl->SetAttribute(el.getAttribute());
      int seid = getHandle().AddElement(newEl);
      m_s2pe.insert({seid, el.getIndex()});

      return *this;
   }

   SubMesh<Context::Serial>& SubMesh<Context::Serial>::add(const BoundaryElement& el)
   {
      assert(&getParent() == &el.getMesh());
      if (isSurface() && !getParent().isSurface())
      {
         assert(getDimension() == getParent().getDimension() - 1);
         assert(getSpaceDimension() == getParent().getSpaceDimension());

         // Add element vertices to the resulting mesh
         mfem::Array<int> pv;
         el.getHandle().GetVertices(pv);

         mfem::Array<int> sv(pv.Size());
         for (int i = 0; i < sv.Size(); i++)
         {
            int pvid = pv[i];
            if (m_s2pv.right.count(pvid) == 0) // Only add vertex if it is not in the map
            {
               sv[i] = getHandle().AddVertex(getParent().getHandle().GetVertex(pvid));
               m_s2pv.insert({sv[i], pvid});
            }
            else // Else get the id of the vertex in the submesh
            {
               sv[i] = m_s2pv.right.at(pvid);
            }
         }

         // Add element with the new vertex ordering
         mfem::Element* pel =
            getHandle().NewElement(el.getHandle().GetGeometryType());
         pel->SetVertices(sv);
         pel->SetAttribute(el.getAttribute());
         int seid = getHandle().AddElement(pel);
         m_s2pe.insert({seid, el.getIndex()});
      }
      else
      {
         assert(getDimension() == getParent().getDimension());
         assert(getSpaceDimension() == getParent().getSpaceDimension());

         // Parent vertices
         mfem::Array<int> pv;
         el.getHandle().GetVertices(pv);

         // Get vertex ordering
         mfem::Array<int> sv(pv.Size());
         for (int i = 0; i < sv.Size(); i++)
            sv[i] = m_s2pv.right.at(pv[i]);

         // Add element with the new vertex ordering
         mfem::Element* newEl =
            getHandle().NewElement(el.getHandle().GetGeometryType());
         newEl->SetVertices(sv);
         newEl->SetAttribute(el.getAttribute());
         int sbid = getHandle().AddBdrElement(newEl);
         m_s2pb.insert({sbid, el.getIndex()});
      }
      return *this;
   }

   const boost::bimap<int, int>& SubMesh<Context::Serial>::getVertexMap() const
   {
      return m_s2pv;
   }

   const boost::bimap<int, int>& SubMesh<Context::Serial>::getElementMap() const
   {
      return m_s2pe;
   }

   const boost::bimap<int, int>& SubMesh<Context::Serial>::getBoundaryElementMap() const
   {
      return m_s2pb;
   }
}
