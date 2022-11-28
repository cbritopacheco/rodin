#ifndef RODIN_MESH_ELEMENT_H
#define RODIN_MESH_ELEMENT_H

#include <set>
#include <mfem.hpp>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
   enum class Type
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

   /**
    * @brief Base class for all geometric elements of the mesh.
    */
   class SimplexBase
   {
      public:
         SimplexBase(const MeshBase& mesh, const mfem::Element* element, int index)
            :  m_mesh(mesh),
               m_element(element),
               m_index(index)
         {}

         SimplexBase(const SimplexBase&) = default;

         SimplexBase(SimplexBase&&) = default;

         virtual ~SimplexBase()
         {
            m_element = nullptr; // This is deallocated by the mfem::Mesh class
         }

         Type getGeometry() const
         {
            return static_cast<Type>(m_element->GetGeometryType());
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
         virtual int getIndex() const
         {
            return m_index;
         }

         virtual std::vector<int> getVertices() const;

         /**
          * @brief Returns the region which the element belongs to.
          */
         virtual Region getRegion() const = 0;

         /**
          * @brief Gets the set of spatially adjacent elements (of the same
          * dimension).
          * @returns Set of integer indices corresponding to the spatially
          * adjacent elements of the same dimension.
          */
         virtual std::set<int> adjacent() const = 0;

         virtual mfem::ElementTransformation& getTransformation() const = 0;

      private:
         const MeshBase& m_mesh;
         const mfem::Element* m_element;
         int m_index;
   };

   bool operator<(const SimplexBase& lhs, const SimplexBase& rhs);

   /**
    * @brief Class for representing elements of the highest dimension in the
    * mesh, i.e. tetrahedra in 3D, triangles in 2D or lines in 1D.
    *
    * This class is designed so that modifications cannot be made to the
    * element. If one wishes to modify the element then one must use
    * ElementView.
    */
   class Element : public SimplexBase
   {
      public:
         Element(const MeshBase& mesh, const mfem::Element* element, int index);

         Element(const Element& other)
            :  SimplexBase(other)
         {
            m_trans.reset(new mfem::IsoparametricTransformation);
            *m_trans = *other.m_trans;
         }

         Element(Element&& other)
            :  SimplexBase(std::move(other)),
               m_trans(std::move(other.m_trans))
         {}

         /**
          * @brief Computes the volume of the element.
          */
         double getVolume() const;

         Region getRegion() const override
         {
            return Region::Domain;
         }

         std::set<int> adjacent() const override;

         virtual mfem::ElementTransformation& getTransformation() const override;

      private:
         std::unique_ptr<mfem::IsoparametricTransformation> m_trans;
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
   class Face : public SimplexBase
   {
      public:
         Face(const MeshBase& mesh, const mfem::Element* element, int index)
            : SimplexBase(mesh, element, index)
         {}

         Face(MeshBase& mesh, mfem::Element* element, int index)
            : SimplexBase(mesh, element, index)
         {}

         Face(const Face& other)
            : SimplexBase(other)
         {}

         Face(Face&& other)
            : SimplexBase(std::move(other))
         {}

         /**
          * @brief Gets the area of the face element.
          */
         double getArea() const;

         std::set<int> elements() const;

         std::set<int> adjacent() const override
         {
            assert(false);
            return {};
         }

         Region getRegion() const override;

         virtual std::set<Element> getElements() const;

         virtual mfem::ElementTransformation& getTransformation() const override;
   };

   class Interface : public Face
   {
      public:
         virtual Region getRegion() const override
         {
            return Region::Interface;
         }
   };

   /**
    * @brief Class for representing boundary elements in the mesh, i.e. face
    * elements which are at the boundary.
    *
    * This class is designed so that modifications can be made to the
    * boundary element. If one wishes to modify the face then one must use
    * BoundaryElementView.
    */
   class Boundary : public Face
   {
      public:
         /**
          *
          * @param[in] index Boundary index
          */
         Boundary(const MeshBase& mesh, const mfem::Element* element, int index);

         Boundary(const Boundary& other)
            :  Face(other),
               m_index(other.m_index)
         {}

         Boundary(Boundary&& other)
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
          */
         int getIndex() const override
         {
            return getBoundaryIndex();
         }

         Region getRegion() const override
         {
            return Region::Boundary;
         }

         mfem::ElementTransformation& getFaceTransformation() const
         {
            return Face::getTransformation();
         }

         mfem::ElementTransformation& getBoundaryTransformation() const
         {
            return *m_trans;
         }

         virtual mfem::ElementTransformation& getTransformation() const override
         {
            return getBoundaryTransformation();
         }

      private:
         int m_index;
         std::unique_ptr<mfem::IsoparametricTransformation> m_trans;
   };

   /**
    * @brief Class for representing boundary elements in the mesh, i.e. face
    * elements which are at the boundary.
    *
    * This class is designed so that modifications cannot be made to the
    * boundary element, while retaining the functionality of the more general
    * FaceView and Boundary classes.
    */
   class BoundaryView
      : public Boundary
   {
      public:
         /**
          * @param[in] index Boundary element index
          */
         BoundaryView(MeshBase& mesh, mfem::Element* element, int index);

         BoundaryView(const BoundaryView& other)
            : Boundary(other),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {}

         BoundaryView(BoundaryView&& other)
            : Boundary(std::move(other)),
              m_mesh(other.m_mesh),
              m_element(other.m_element)
         {
            other.m_element = nullptr;
         }

         /**
          * @brief Sets the attribute of the boundary element.
          */
         BoundaryView& setAttribute(int attr);

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
    * @brief Represents a spatial point which belongs to some element of a mesh.
    *
    * A Point represents the physical coordinates (as opposed to reference
    * coordinates) of a point on the mesh.
    * This class differs from a Vertex in the sense that a Vertex is a node of
    * some element of a Mesh. In contrast, a Point represents any
    * coordinates contained in the Mesh.
    */
   class Point
   {
      public:
         Point(const mfem::IntegrationPoint& ip)
            : m_ip(ip)
         {}

         /**
          * @brief Constructs the Point object from reference coordinates.
          * @param[in] trans Element transformation
          * @param[in] ip Reference coordinates
          */
         Point(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
            :  m_trans(trans),
               m_ip(ip)
         {
            assert(&trans.GetIntPoint() == &ip);
            m_trans->get().Transform(ip, m_physical);
            assert(m_physical.Size() == trans.mesh->SpaceDimension());
         }

         Point(const Point&) = default;

         Point(Point&&) = default;

         bool transformed() const
         {
            return m_trans.has_value();
         }

         Point& transform(mfem::ElementTransformation& trans)
         {
            trans.SetIntPoint(&m_ip.get());
            m_trans.emplace(trans);
            m_trans->get().Transform(m_ip, m_physical);
            return *this;
         }

         /**
          * @brief Gets the space dimension of the physical coordinates.
          * @returns Dimension of the physical coordinates.
          */
         int getDimension() const
         {
            assert(m_trans);
            return m_trans->get().mesh->SpaceDimension();
         }

         /**
          * @brief Gets the i-th physical coordinate.
          * @returns Physical i-th coordinate.
          */
         double operator()(int i) const
         {
            return m_physical(i);
         }

         /**
          * @brief Gets the @f$ x @f$ physical coordinate.
          * @returns Physical @f$ x @f$-coordinate.
          */
         double x() const
         {
            return m_physical(0);
         }

         /**
          * @brief Gets the @f$ y @f$ physical coordinate.
          * @returns Physical @f$ y @f$-coordinate.
          */
         double y() const
         {
            return m_physical(1);
         }

         /**
          * @brief Gets the @f$ z @f$ physical coordinate.
          * @returns Physical @f$ z @f$-coordinate.
          */
         double z() const
         {
            return m_physical(2);
         }

         /**
          * @brief Lexicographical comparison.
          */
         bool operator<(const Point& rhs) const
         {
            assert(getDimension() == rhs.getDimension());
            for (int i = 0; i < m_physical.Size() - 1; i++)
            {
               if (m_physical(i) < rhs.m_physical(i))
                  return true;
               if (rhs.m_physical(i) > m_physical(i))
                  return false;
            }
            return (m_physical(m_physical.Size() - 1) < rhs.m_physical(rhs.m_physical.Size() - 1));
         }

         mfem::ElementTransformation& getElementTransformation() const
         {
            assert(m_trans);
            return m_trans.value().get();
         }

         const mfem::IntegrationPoint& getIntegrationPoint() const
         {
            return m_ip.get();
         }

      private:
         mfem::Vector m_physical;
         std::optional<std::reference_wrapper<mfem::ElementTransformation>> m_trans;
         std::reference_wrapper<const mfem::IntegrationPoint> m_ip;
   };
}

#endif
