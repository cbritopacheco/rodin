#include "Element.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
   bool operator<(const Simplex& lhs, const Simplex& rhs)
   {
      return lhs.getIndex() < rhs.getIndex();
   }

   // ---- Simplex -----------------------------------------------------------
   Index Simplex::getIndex() const
   {
      return m_data.index;
   }

   Type Simplex::getGeometry() const
   {
      return static_cast<Type>(m_data.element->GetType());
   }

   Attribute Simplex::getAttribute() const
   {
      return m_data.element->GetAttribute();
   }

   const MeshBase& Simplex::getMesh() const
   {
      return m_data.mesh.get();
   }

   mfem::ElementTransformation& Simplex::getTransformation() const
   {
      return *m_data.trans;
   }

   std::vector<Geometry::Point> Simplex::getIntegrationRule(int order) const
   {
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(getTransformation().GetGeometryType(), order);
      std::vector<Geometry::Point> res;
      res.reserve(ir->GetNPoints());
      for (int i = 0; i < ir->GetNPoints(); i++)
         res.emplace_back(*this, ir->IntPoint(i));
      return res;
   }

   double Simplex::getVolume() const
   {
      mfem::ElementTransformation& trans = getTransformation();
      const mfem::IntegrationRule& ir =
         mfem::IntRules.Get(static_cast<int>(getGeometry()), trans.OrderJ());
      double volume = 0.0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(j);
         trans.SetIntPoint(&ip);
         volume += ip.weight * trans.Weight();
      }
      return volume;
   }

   SimplexVertexIterator Simplex::getVertices() const
   {
      // mfem::Array<int> vs;
      // m_data.element->GetVertices(vs);
      // return std::vector<int>(vs.begin(), vs.end());
      assert(false);
   }

   // ---- Element -----------------------------------------------------------
   AdjacentElementIterator Element::getAdjacent() const
   {
      assert(false);
   }

   // ---- Face --------------------------------------------------------------
   Face::Face(Data data)
      : Simplex(std::move(data))
   {}

   FaceElementIterator Face::getIncident() const
   {
      assert(false);
   }

   Region Face::getRegion() const
   {
      if (getMesh().isInterface(getIndex()))
         return Region::Interface;
      else
         return Region::Boundary;
   }

   // std::set<int> Face::elements() const
   // {
   //    int e1 = -1, e2 = -1;
   //    getMesh().getHandle().GetFaceElements(Face::getIndex(), &e1, &e2);
   //    if (e1 >= 0 && e2 >= 0)
   //       return {e1, e2};
   //    else if (e1 >= 0 && e2 < 0)
   //       return {e1};
   //    else if (e1 < 0 && e2 >= 0)
   //       return {e2};
   //    else
   //       return {};
   // }

   // std::set<Element> Face::getElements() const
   // {
   //    int e1 = -1, e2 = -1;
   //    getMesh().getHandle().GetFaceElements(Face::getIndex(), &e1, &e2);
   //    if (e1 >= 0 && e2 >= 0)
   //       return { Element(getMesh(), getMesh().getHandle().GetElement(e1), e1),
   //                Element(getMesh(), getMesh().getHandle().GetElement(e2), e2) };
   //    else if (e1 >= 0 && e2 < 0)
   //       return { Element(getMesh(), getMesh().getHandle().GetElement(e1), e1) };
   //    else if (e1 < 0 && e2 >= 0)
   //       return { Element(getMesh(), getMesh().getHandle().GetElement(e2), e2) };
   //    else
   //       return {};
   // }

   // ---- Point -------------------------------------------------------------
   Point::Point(const Simplex& element, const mfem::IntegrationPoint& ip)
      : m_element(element), m_ip(ip)
   {
      m_element.get().getTransformation().SetIntPoint(&m_ip);
      m_element.get().getTransformation().Transform(m_ip, m_physical);
      assert(static_cast<size_t>(m_physical.Size()) == element.getMesh().getSpaceDimension());
   }

   Point::Point(const Simplex& element, mfem::IntegrationPoint&& ip)
      : m_element(element), m_ip(std::move(ip))
   {}

   int Point::getDimension(Coordinates coords) const
   {
      switch (coords)
      {
         case Coordinates::Physical:
         {
            return m_element.get().getMesh().getSpaceDimension();
         }
         case Coordinates::Reference:
         {
            return m_element.get().getMesh().getDimension();
         }
      }
   }
}
