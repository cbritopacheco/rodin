#include "SubMesh.h"

namespace Rodin
{
   SubMesh<Traits::Serial>::SubMesh(const SubMesh& other)
      :  Mesh(other),
         m_parent(other.m_parent),
         m_s2pv(other.m_s2pv)
   {}

   const boost::bimap<int, int>& SubMesh<Traits::Serial>::getVertexMap() const
   {
      assert(m_s2pv);
      return *m_s2pv;
   }

   SubMesh<Traits::Serial>& SubMesh<Traits::Serial>::setVertexMap(boost::bimap<int, int> s2pv)
   {
      m_s2pv = s2pv;
      return *this;
   }

   SubMesh<Traits::Serial>& SubMesh<Traits::Serial>::setParent(const MeshBase& parent)
   {
      m_parent = std::cref(parent);
      return *this;
   }

   const MeshBase& SubMesh<Traits::Serial>::getParent() const
   {
      assert(m_parent);
      return m_parent->get();
   }
}
