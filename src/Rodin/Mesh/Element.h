#ifndef RODIN_MESH_ELEMENT_H
#define RODIN_MESH_ELEMENT_H

#include <mfem.hpp>
#include "ForwardDecls.h"

namespace Rodin
{
   class ElementBase
   {
      public:
         enum Geometry
         {
            Invalid = mfem::Geometry::INVALID,
            Point = mfem::Geometry::POINT,
            Segment = mfem::Geometry::SEGMENT,
            Triangle = mfem::Geometry::TRIANGLE,
            Square = mfem::Geometry::SQUARE,
            Tetrahedron = mfem::Geometry::TETRAHEDRON,
            Cube = mfem::Geometry::CUBE,
            Prism = mfem::Geometry::PRISM,
            Pyramid = mfem::Geometry::PYRAMID
         };

         ElementBase(MeshBase& mesh, mfem::Element* element, int index);

         virtual ~ElementBase()
         {
            m_element = nullptr; // This is deallocated by the mfem::Mesh class
         }

         Geometry getGeometry() const
         {
            return static_cast<Geometry>(m_element->GetGeometryType());
         }

         double getAttribute() const;

         int getIndex() const
         {
            return m_index;
         }

         MeshBase& getMesh()
         {
            return m_mesh;
         }

         const MeshBase& getMesh() const
         {
            return m_mesh;
         }

         mfem::Element& getHandle()
         {
            return *m_element;
         }

         const mfem::Element& getHandle() const
         {
            return *m_element;
         }

      private:
         MeshBase& m_mesh;
         mfem::Element* m_element;
         int m_index;
   };

   class Element : public ElementBase
   {
      public:
         Element(MeshBase& mesh, mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         double getVolume() const;
   };

   class BoundaryElement : public ElementBase
   {
      public:
         BoundaryElement(MeshBase& mesh, mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         double getArea() const;
   };
}

#endif
