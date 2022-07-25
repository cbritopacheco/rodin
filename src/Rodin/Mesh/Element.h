#ifndef RODIN_MESH_ELEMENT_H
#define RODIN_MESH_ELEMENT_H

#include <set>
#include <mfem.hpp>

#include "ForwardDecls.h"

namespace Rodin
{
   /**
    * @brief Base class for all geometric elements of the mesh.
    */
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

         ElementBase(const ElementBase&) = default;

         ElementBase(ElementBase&&) = default;

         virtual ~ElementBase()
         {
            m_element = nullptr; // This is deallocated by the mfem::Mesh class
         }

         Geometry getGeometry() const
         {
            return static_cast<Geometry>(m_element->GetGeometryType());
         }

         /**
          * @brief Gets the attribute of the element.
          */
         int getAttribute() const;

         /**
          * @brief Gets the index of the element in the mesh.
          */
         int getIndex() const
         {
            return m_index;
         }

         /**
          * @brief Gets the associated mesh to the element.
          */
         const MeshBase& getMesh() const
         {
            return m_mesh;
         }

         const mfem::Element& getHandle() const
         {
            return *m_element;
         }

         /**
          * @brief Gets the set of spatially adjacent elements (of the same
          * dimension).
          * @returns Set of integer indices corresponding to the spatially
          * adjacent elements of the same dimension.
          */
         virtual std::set<int> adjacent() const = 0;

      private:
         const MeshBase& m_mesh;
         const mfem::Element* m_element;
         int m_index;
   };

   /**
    * @brief Class for representing elements of the highest dimension in the
    * mesh, i.e. tetrahedra in 3D, triangles in 2D or lines in 1D.
    *
    * This class is designed so that modifications cannot be made to the
    * element. If one wishes to modify the element then one must use
    * ElementView.
    */
   class Element : public ElementBase
   {
      public:
         Element(const MeshBase& mesh, const mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         Element(const Element& other)
            : ElementBase(other)
         {}

         Element(Element&& other)
            : ElementBase(std::move(other))
         {}

         /**
          * @brief Computes the volume of the element.
          */
         double getVolume() const;

         std::set<int> adjacent() const override;
   };

   /**
    * @brief Class for representing elements of the highest dimension in the
    * mesh, i.e. tetrahedra in 3D, triangles in 2D or lines in 1D.
    *
    * This class is designed so that modifications can be made to the element,
    * while retaining the functionality of the more general Element class.
    */
   class ElementView : public Element
   {
      public:
         ElementView(MeshBase& mesh, mfem::Element* element, int index)
            : Element(mesh, element, index),
              m_mesh(mesh),
              m_element(element)
         {}

         ElementView(const ElementView& other)
            : Element(other),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {}

         ElementView(ElementView&& other)
            : Element(std::move(other)),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {
            other.m_element = nullptr;
         }

         /**
          * @brief Sets the attribute of the element.
          */
         ElementView& setAttribute(int attr);

         MeshBase& getMesh()
         {
            return m_mesh;
         }

         mfem::Element& getHandle()
         {
            return *m_element;
         }

      private:
         MeshBase& m_mesh;
         mfem::Element* m_element;

   };

   /**
    * @brief Class for representing elements of codimension 1 in the
    * mesh, i.e. triangles in 3D or lines in 2D.
    *
    * This class is designed so that modifications cannot be made to the
    * face. If one wishes to modify the face then one must use
    * FaceView.
    */
   class Face : public ElementBase
   {
      public:
         Face(const MeshBase& mesh, const mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         Face(MeshBase& mesh, mfem::Element* element, int index)
            : ElementBase(mesh, element, index)
         {}

         Face(const Face& other)
            : ElementBase(other)
         {}

         Face(Face&& other)
            : ElementBase(std::move(other))
         {}

         std::set<int> adjacent() const override
         {
            assert(false);
         }

         /**
          * @brief Gets the area of the face element.
          */
         double getArea() const;
   };

   /**
    * @brief Class for representing elements of codimension 1 in the
    * mesh, i.e. triangles in 3D or lines in 2D.
    *
    * This class is designed so that modifications can be made to the face,
    * while retaining the functionality of the more general face class.
    */
   class FaceView : public Face
   {
      public:
         FaceView(MeshBase& mesh, mfem::Element* element, int index)
            : Face(mesh, element, index),
              m_mesh(mesh),
              m_element(element)
         {}

         FaceView(const FaceView& other)
            : Face(other),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {}

         FaceView(FaceView&& other)
            : Face(std::move(other)),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {
            m_element = nullptr;
         }

         MeshBase& getMesh()
         {
            return m_mesh;
         }

         mfem::Element& getHandle()
         {
            return *m_element;
         }

      private:
         MeshBase& m_mesh;
         mfem::Element* m_element;
   };

   /**
    * @brief Class for representing boundary elements in the mesh, i.e. face
    * elements which are at the boundary.
    *
    * This class is designed so that modifications can be made to the
    * boundary element. If one wishes to modify the face then one must use
    * BoundaryElementView.
    */
   class BoundaryElement : public Face
   {
      public:
         BoundaryElement(const MeshBase& mesh, const mfem::Element* element, int index)
            : Face(mesh, element, index)
         {}

         BoundaryElement(const BoundaryElement& other)
            :  Face(other)
         {}

         BoundaryElement(BoundaryElement&& other)
            :  Face(std::move(other))
         {}
   };

   /**
    * @brief Class for representing boundary elements in the mesh, i.e. face
    * elements which are at the boundary.
    *
    * This class is designed so that modifications cannot be made to the
    * boundary element, while retaining the functionality of the more general
    * FaceView and BoundaryElement classes.
    */
   class BoundaryElementView : public FaceView, public BoundaryElement
   {
      public:
         BoundaryElementView(MeshBase& mesh, mfem::Element* element, int index)
            : FaceView(mesh, element, index),
              BoundaryElement(mesh, element, index)
         {}

         BoundaryElementView(const BoundaryElementView& other)
            :  FaceView(other),
               BoundaryElement(other)
         {}

         BoundaryElementView(BoundaryElementView&& other)
            :  FaceView(std::move(other)),
               BoundaryElement(std::move(other))
         {}

         int getAttribute() const
         {
            return BoundaryElement::getAttribute();
         }

         /**
          * @brief Sets the attribute of the boundary element.
          */
         BoundaryElementView& setAttribute(int attr);
   };
}

#endif
