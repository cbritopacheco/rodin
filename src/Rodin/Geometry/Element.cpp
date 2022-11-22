#include "Element.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
   bool operator<(const SimplexBase& lhs, const SimplexBase& rhs)
   {
      return lhs.getIndex() < rhs.getIndex();
   }

   // ---- ElementBase -------------------------------------------------------
   std::vector<int> SimplexBase::getVertices() const
   {
      mfem::Array<int> vs;
      m_element->GetVertices(vs);
      return std::vector<int>(vs.begin(), vs.end());
   }

   int SimplexBase::getAttribute() const
   {
      return m_element->GetAttribute();
   }

   // ---- Element -----------------------------------------------------------
   std::set<int> Element::adjacent() const
   {
      std::set<int> res;
      // This call does not actually modify the Element object. Only the Mesh
      // object. Hence we const_cast to access the ElementToElementTable table
      // method.
      const auto& adj = const_cast<MeshBase&>(getMesh()).getHandle().ElementToElementTable();
      const int* elements = adj.GetRow(getIndex());
      for (int i = 0; i < adj.RowSize(getIndex()); i++)
         res.insert(elements[i]);
      return res;
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

   mfem::ElementTransformation& Element::getTransformation() const
   {
      return *const_cast<mfem::Mesh&>(getMesh().getHandle()).GetElementTransformation(getIndex());
   }

   // ---- ElementView -------------------------------------------------------
   ElementView& ElementView::setAttribute(int attr)
   {
      getMesh().getHandle().SetAttribute(getIndex(), attr);
      return *this;
   }

   // ---- Face --------------------------------------------------------------
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

   std::set<int> Face::elements() const
   {
      int e1 = -1, e2 = -1;
      getMesh().getHandle().GetFaceElements(Face::getIndex(), &e1, &e2);
      if (e1 >= 0 && e2 >= 0)
         return {e1, e2};
      else if (e1 >= 0 && e2 < 0)
         return {e1};
      else if (e1 < 0 && e2 >= 0)
         return {e2};
      else
         return {};
   }

   std::set<Element> Face::getElements() const
   {
      int e1 = -1, e2 = -1;
      getMesh().getHandle().GetFaceElements(Face::getIndex(), &e1, &e2);
      if (e1 >= 0 && e2 >= 0)
         return { Element(getMesh(), getMesh().getHandle().GetElement(e1), e1),
                  Element(getMesh(), getMesh().getHandle().GetElement(e2), e2) };
      else if (e1 >= 0 && e2 < 0)
         return { Element(getMesh(), getMesh().getHandle().GetElement(e1), e1) };
      else if (e1 < 0 && e2 >= 0)
         return { Element(getMesh(), getMesh().getHandle().GetElement(e2), e2) };
      else
         return {};
   }

   Region Face::getRegion() const
   {
      if (getMesh().getHandle().FaceIsInterior(getIndex()))
         return Region::Interface;
      else
         return Region::Boundary;
   }

   mfem::ElementTransformation& Face::getTransformation() const
   {
      return *const_cast<mfem::Mesh&>(getMesh().getHandle()).GetFaceTransformation(Face::getIndex());
   }

   mfem::FaceElementTransformations& Face::getTransformations() const
   {
      return *const_cast<mfem::Mesh&>(getMesh().getHandle()).GetFaceElementTransformations(Face::getIndex());
   }

   // ---- BoundaryElement ---------------------------------------------------
   Boundary::Boundary(
         const MeshBase& mesh, const mfem::Element* element, int index)
      : Face(mesh, element, mesh.getHandle().GetBdrFace(index)),
        m_index(index)
   {}

   mfem::ElementTransformation& Boundary::getTransformation() const
   {
      return *const_cast<mfem::Mesh&>(getMesh().getHandle()).GetBdrElementTransformation(getBoundaryIndex());
   }

   // ---- BoundaryElementView -----------------------------------------------
   BoundaryView::BoundaryView(MeshBase& mesh, mfem::Element* element, int index)
      : Boundary(mesh, element, index),
        m_mesh(mesh),
        m_element(element)
   {}

   BoundaryView& BoundaryView::setAttribute(int attr)
   {
      getMesh().getHandle().SetBdrAttribute(getIndex(), attr);
      return *this;
   }

   // ---- Vertex ------------------------------------------------------------
}
