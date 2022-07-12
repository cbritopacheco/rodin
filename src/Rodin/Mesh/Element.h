#ifndef RODIN_MESH_ELEMENT_H
#define RODIN_MESH_ELEMENT_H

#include <set>
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

         ElementBase(const MeshBase& mesh, const mfem::Element* element, int index)
            :  m_mesh(mesh),
               m_element(element),
               m_index(index)
         {}

         virtual ~ElementBase()
         {
            m_element = nullptr; // This is deallocated by the mfem::Mesh class
         }

         Geometry getGeometry() const
         {
            return static_cast<Geometry>(m_element->GetGeometryType());
         }

         int getAttribute() const;

         int getIndex() const
         {
            return m_index;
         }

         const MeshBase& getMesh() const
         {
            return m_mesh;
         }

         const mfem::Element& getHandle() const
         {
            return *m_element;
         }

         virtual std::set<int> adjacent() const = 0;

      private:
         const MeshBase& m_mesh;
         const mfem::Element* m_element;
         int m_index;
   };

   class Element : public ElementBase
   {
      public:
         Element(const MeshBase& mesh, const mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         double getVolume() const;

         std::set<int> adjacent() const override;
   };

   class ElementView : public Element
   {
      public:
         ElementView(MeshBase& mesh, mfem::Element* element, int index)
            : Element(mesh, element, index),
              m_mesh(mesh),
              m_element(element)
         {}

         MeshBase& getMesh()
         {
            return m_mesh;
         }

         mfem::Element* getHandle()
         {
            return m_element;
         }

      private:
         MeshBase& m_mesh;
         mfem::Element* m_element;

   };

   class Face : public ElementBase
   {
      public:
         Face(const MeshBase& mesh, const mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         std::set<int> adjacent() const override
         {
            assert(false);
         }

         double getArea() const;
   };

   class FaceView : public Face
   {
      public:
         FaceView(MeshBase& mesh, mfem::Element* element, int index)
            : Face(mesh, element, index),
              m_mesh(mesh),
              m_element(element)
         {}

         MeshBase& getMesh()
         {
            return m_mesh;
         }

         mfem::Element* getHandle()
         {
            return m_element;
         }

      private:
         MeshBase& m_mesh;
         mfem::Element* m_element;
   };
}

#endif
