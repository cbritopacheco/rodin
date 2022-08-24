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
          * @brief Gets the index of the element in the mesh.
          */
         int getIndex() const
         {
            return m_index;
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
    * face.
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

         std::set<int> elements() const;

         std::set<int> adjacent() const override
         {
            assert(false);
            return {};
         }

         /**
          * @brief Gets the area of the face element.
          */
         double getArea() const;
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
         /**
          *
          * @param[in] index Boundary index
          */
         BoundaryElement(
               const MeshBase& mesh, const mfem::Element* element, int index);

         BoundaryElement(const BoundaryElement& other)
            :  Face(other),
               m_index(other.m_index)
         {}

         BoundaryElement(BoundaryElement&& other)
            :  Face(std::move(other)),
               m_index(other.m_index)
         {
            other.m_index = 0;
         }

         /**
          * @brief Gets the index of the face element.
          *
          * @warning The index returned by this function need not coincide with
          * the index returned by getBoundaryIndex().
          */
         int getFaceIndex() const
         {
            return Face::getIndex();
         }

         /**
          * @brief Gets the index of the element as a member of the boundary.
          *
          * @warning The index returned by this function need not coincide with
          * the index returned by getFaceIndex().
          */
         int getBoundaryIndex() const
         {
            return m_index;
         }

         /**
          * @brief Gets the index of the element as a member of the boundary.
          *
          * Same as getBoundaryIndex().
          *
          * @warning The index returned by this function need not coincide with
          * the index returned by getFaceIndex().
          */
         int getIndex() const
         {
            return getBoundaryIndex();
         }

      private:
         int m_index;
   };

   /**
    * @brief Class for representing boundary elements in the mesh, i.e. face
    * elements which are at the boundary.
    *
    * This class is designed so that modifications cannot be made to the
    * boundary element, while retaining the functionality of the more general
    * FaceView and BoundaryElement classes.
    */
   class BoundaryElementView
      : public BoundaryElement
   {
      public:
         /**
          * @param[in] index Boundary element index
          */
         BoundaryElementView(MeshBase& mesh, mfem::Element* element, int index);

         BoundaryElementView(const BoundaryElementView& other)
            : BoundaryElement(other),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {}

         BoundaryElementView(BoundaryElementView&& other)
            : BoundaryElement(std::move(other)),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {
            other.m_element = nullptr;
         }

         /**
          * @brief Sets the attribute of the boundary element.
          */
         BoundaryElementView& setAttribute(int attr);

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

   class Vertex
   {
      public:
         Vertex(mfem::Vector&& v)
            : m_v(std::move(v)),
              m_trans(nullptr),
              m_ip(nullptr)
         {}

         Vertex(const mfem::Vector& v)
            : m_v(v),
              m_trans(nullptr),
              m_ip(nullptr)
         {}

         Vertex(const Vertex& other) = default;

         Vertex(Vertex&& other)
            : m_v(std::move(other.m_v)),
              m_trans(other.m_trans),
              m_ip(other.m_ip)
         {
            m_trans = nullptr;
            m_ip = nullptr;
         }

         int getDimension() const
         {
            return m_v.Size();
         }

         double& operator()(int i)
         {
            return m_v(i);
         }

         const double& operator()(int i) const
         {
            return m_v(i);
         }

         double& x()
         {
            return m_v(0);
         }

         double& y()
         {
            return m_v(1);
         }

         double& z()
         {
            return m_v(2);
         }

         const double& x() const
         {
            return m_v(0);
         }

         const double& y() const
         {
            return m_v(1);
         }

         const double& z() const
         {
            return m_v(2);
         }

         bool operator<(const Vertex& rhs) const
         {
            assert(getDimension() == rhs.getDimension());
            bool r = true;
            for (int i = 0; i < m_v.Size(); i++)
               r = r && m_v(i) < rhs.m_v(i);
            return r;
         }

         Vertex& setElementTransformation(mfem::ElementTransformation* trans)
         {
            m_trans = trans;
            return *this;
         }

         Vertex& setIntegrationPoint(const mfem::IntegrationPoint* ip)
         {
            m_ip = ip;
            return *this;
         }

         mfem::ElementTransformation* getElementTransformation() const
         {
            return m_trans;
         }

         const mfem::IntegrationPoint* getIntegrationPoint() const
         {
            return m_ip;
         }

      private:
         mfem::Vector m_v;
         mfem::ElementTransformation* m_trans;
         const mfem::IntegrationPoint* m_ip;
   };
}

#endif
