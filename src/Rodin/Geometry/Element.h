/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_ELEMENT_H
#define RODIN_MESH_ELEMENT_H

#include <set>
#include <array>
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
   class Simplex
   {
      public:
         struct Data
         {
            std::unique_ptr<mfem::Element> element;
            std::unique_ptr<mfem::ElementTransformation> trans;
         };

         Simplex(size_t dimension, Index index, const MeshBase& mesh, Data data)
            :  m_dimension(dimension), m_index(index), m_mesh(mesh),
               m_data(std::move(data))
         {}

         Simplex(Simplex&&) = default;

         virtual ~Simplex() = default;

         size_t getDimension() const
         {
            return m_dimension;
         }

         /**
          * @brief Gets the index of the simplex in the mesh.
          */
         Index getIndex() const
         {
            return m_index;
         }

         /**
          * @brief Gets the associated mesh to the simplex.
          */
         const MeshBase& getMesh() const
         {
            return m_mesh.get();
         }

         Type getGeometry() const;

         /**
          * @brief Gets the attribute of the simplex.
          */
         Attribute getAttribute() const;

         double getVolume() const;

         // VertexIterator getVertices() const;

         mfem::ElementTransformation& getTransformation() const;

         virtual std::vector<Geometry::Point> getIntegrationRule(int order) const;

         const mfem::Element& getHandle() const
         {
            assert(m_data.element);
            return *m_data.element;
         }

      private:
         const size_t m_dimension;
         const Index m_index;
         std::reference_wrapper<const MeshBase> m_mesh;
         Data m_data;
   };

   bool operator<(const Simplex& lhs, const Simplex& rhs);

   /**
    * @brief Class for representing elements of the highest dimension in the
    * mesh, i.e. tetrahedra in 3D, triangles in 2D or lines in 1D.
    *
    * This class is designed so that modifications cannot be made to the
    * element. If one wishes to modify the element then one must use
    * ElementView.
    */
   class Element : public Simplex
   {
      public:
         Element(Index index, const MeshBase& mesh, Data data);

         Element(Element&& other)
            :  Simplex(std::move(other))
         {}

         ElementIterator getAdjacent() const;
   };

   /**
    * @brief Class for representing elements of codimension 1 in the
    * mesh, i.e. triangles in 3D or lines in 2D.
    *
    * This class is designed so that modifications cannot be made to the
    * face.
    */
   class Face : public Simplex
   {
      public:
         Face(Index index, const MeshBase& mesh, Data data);

         Face(Face&& other)
            : Simplex(std::move(other))
         {}

         bool isBoundary() const;

         bool isInterface() const;

         FaceIterator getAdjacent() const;

         ElementIterator getIncident() const;

         mfem::FaceElementTransformations& getFaceTransformations() const
         {
            assert(m_localTrans);
            return *m_localTrans;
         }

      private:
         std::unique_ptr<mfem::FaceElementTransformations> m_localTrans;
   };

   class Interface final : public Face
   {
      public:
         Interface(Index index, const MeshBase& mesh, Data data)
            : Face(index, mesh, std::move(data))
         {}

         Interface(Interface&& other)
            : Face(std::move(other))
         {}
   };

   /**
    * @brief Class for representing boundary elements in the mesh, i.e. face
    * elements which are at the boundary.
    *
    * This class is designed so that modifications can be made to the
    * boundary element. If one wishes to modify the face then one must use
    * BoundaryElementView.
    */
   class Boundary final : public Face
   {
      public:
         Boundary(Index index, const MeshBase& mesh, Data data)
            : Face(index, mesh, std::move(data))
         {}

         Boundary(Boundary&& other)
            :  Face(std::move(other))
         {}
   };

   class Vertex : public Simplex
   {
      public:
         double x() const
         {
            return operator()(0);
         }

         double y() const
         {
            return operator()(1);
         }

         double z() const
         {
            return operator()(2);
         }

         virtual double operator()(size_t i) const;
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
         enum class Coordinates
         {
            Reference,
            Physical
         };

         /**
          * @brief Constructs the Point object from reference coordinates.
          * @param[in] simplex Simplex to which point belongs to
          * @param[in] ip Reference coordinates
          */
         Point(const Simplex& simplex, const mfem::IntegrationPoint& ip);

         Point(const Simplex& simplex, mfem::IntegrationPoint&& ip);

         Point(const Point&) = default;

         Point(Point&&) = default;

         /**
          * @brief Gets the space dimension of the physical coordinates.
          * @returns Dimension of the physical coordinates.
          */
         int getDimension(Coordinates coords = Coordinates::Physical) const;

         /**
          * @brief Gets the i-th physical coordinate.
          * @returns Physical i-th coordinate.
          */
         double operator()(int i, Coordinates coords = Coordinates::Physical) const
         {
            switch (coords)
            {
               case Coordinates::Physical:
               {
                  assert(m_physical.Size() > i);
                  return m_physical(i);
               }
               case Coordinates::Reference:
               {
                  double p[3];
                  m_ip.Get(p, getDimension(coords));
                  return p[i];
               }
            }
         }

         /**
          * @brief Gets the @f$ x @f$ physical coordinate.
          * @returns Physical @f$ x @f$-coordinate.
          */
         double x(Coordinates coords = Coordinates::Physical) const
         {
            switch (coords)
            {
               case Coordinates::Physical:
               {
                  assert(m_physical.Size() > 0);
                  return m_physical(0);
               }
               case Coordinates::Reference:
               {
                  return m_ip.x;
               }
            }
         }

         /**
          * @brief Gets the @f$ y @f$ physical coordinate.
          * @returns Physical @f$ y @f$-coordinate.
          */
         double y(Coordinates coords = Coordinates::Physical) const
         {
            switch (coords)
            {
               case Coordinates::Physical:
               {
                  assert(m_physical.Size() > 1);
                  return m_physical(1);
               }
               case Coordinates::Reference:
               {
                  return m_ip.y;
               }
            }
         }

         /**
          * @brief Gets the @f$ z @f$ physical coordinate.
          * @returns Physical @f$ z @f$-coordinate.
          */
         double z(Coordinates coords = Coordinates::Physical) const
         {
            switch (coords)
            {
               case Coordinates::Physical:
               {
                  assert(m_physical.Size() > 2);
                  return m_physical(2);
               }
               case Coordinates::Reference:
               {
                  return m_ip.z;
               }
            }
         }

         double w() const
         {
            return m_ip.weight;
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

         const Simplex& getSimplex() const
         {
            return m_element;
         }

      private:
         mfem::Vector m_physical;
         std::reference_wrapper<const Simplex> m_element;
         mfem::IntegrationPoint m_ip;
   };
}

#endif
