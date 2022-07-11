#include "Element.h"
#include "Mesh.h"

namespace Rodin
{
   ElementBase::ElementBase(MeshBase& mesh, mfem::Element* element, int index)
      : m_mesh(mesh), m_element(element), m_index(index)
   {}

   double ElementBase::getAttribute() const
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

   double BoundaryElement::getArea() const
   {
      mfem::ElementTransformation *et = const_cast<MeshBase&>(getMesh()
            ).getHandle().GetBdrElementTransformation(getIndex());
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
