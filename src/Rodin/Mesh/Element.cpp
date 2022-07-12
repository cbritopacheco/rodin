#include "Element.h"
#include "Mesh.h"

namespace Rodin
{
   std::set<int> Element::adjacent() const
   {
      std::set<int> res;
      // This call does not actually modify the Element object. Only the Mesh
      // object. Hence we const_cast to access the ElementToElementTable table
      // method.
      auto& adj = const_cast<MeshBase&>(getMesh()).getHandle().ElementToElementTable();
      const int* elements = adj.GetRow(getIndex());
      for (int i = 0; i < adj.RowSize(getIndex()); i++)
         res.insert(elements[i]);
      return res;
   }

   int ElementBase::getAttribute() const
   {
      return m_element->GetAttribute();
   }

   double Element::getVolume() const
   {
      mfem::ElementTransformation *et = const_cast<MeshBase&>(getMesh()
            ).getHandle().GetElementTransformation(getIndex());
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(
            getHandle().GetGeometryType(), et->OrderJ());
      double volume = 0.0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(j);
         et->SetIntPoint(&ip);
         volume += ip.weight * et->Weight();
      }
      return volume;
   }

   double Face::getArea() const
   {
      mfem::ElementTransformation *et = const_cast<MeshBase&>(getMesh()
            ).getHandle().GetFaceTransformation(getIndex());
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(
            getHandle().GetGeometryType(), et->OrderJ());
      double area = 0.0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(j);
         et->SetIntPoint(&ip);
         area += ip.weight * et->Weight();
      }
      return area;
   }
}
