#include "SubMesh.h"

#include "Element.h"

namespace Rodin
{
   SubMesh<Traits::Serial>::SubMesh(const MeshBase& parent)
      : m_parent(parent)
   {}

   SubMesh<Traits::Serial>::SubMesh(const SubMesh& other)
      :  Mesh(other),
         m_parent(other.m_parent),
         m_s2pv(other.m_s2pv)
   {}

   const boost::bimap<int, int>& SubMesh<Traits::Serial>::getVertexMap() const
   {
      return m_s2pv;
   }

   SubMesh<Traits::Serial>& SubMesh<Traits::Serial>::setVertexMap(boost::bimap<int, int> s2pv)
   {
      m_s2pv = s2pv;
      return *this;
   }

   const MeshBase& SubMesh<Traits::Serial>::getParent() const
   {
      return m_parent;
   }

   SubMesh<Traits::Serial>& SubMesh<Traits::Serial>::add(const ElementBase& el)
   {
      assert(&getParent() == &el.getMesh());

      // Add element vertices to the resulting mesh
      mfem::Array<int> pv;
      el.getHandle().GetVertices(pv);

      mfem::Array<int> sv(pv.Size());
      for (int i = 0; i < sv.Size(); i++)
      {
         int pvid = pv[i];
         if (m_s2pv.right.count(pvid) == 0) // Only add vertex if it is not in the map
            sv[i] = getHandle().AddVertex(getParent().getHandle().GetVertex(pvid));
      }

      // Add element with the new vertex ordering
      mfem::Element* pel =
         getHandle().NewElement(el.getHandle().GetGeometryType());
      pel->SetVertices(sv);
      pel->SetAttribute(el.getAttribute());
      int seid = getHandle().AddElement(pel);
      m_s2pe.insert({seid, el.getIndex()});

      return *this;
   }
}
