/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEX_H
#define RODIN_GEOMETRY_SIMPLEX_H

#include <set>
#include <array>
#include <mfem.hpp>

#include "Rodin/Math/Vector.h"

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
      enum class Property
      {
        Attribute
      };

      Simplex(
          size_t dimension,
          Index index,
          const MeshBase& mesh,
          const std::vector<Index>& vertices,
          Attribute attr = RODIN_DEFAULT_SIMPLEX_ATTRIBUTE);

      Simplex(const Simplex&) = delete;

      Simplex(Simplex&&) = default;

      virtual ~Simplex() = default;

      /**
       * @brief Gets the index of the simplex in the mesh.
       */
      Index getIndex() const
      {
        return m_index;
      }


      Type getGeometry() const
      {
        return m_type;
      }

      double getVolume() const;

      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the attribute of the simplex.
       */
      Attribute getAttribute() const
      {
        return m_attr;
      }

      /**
       * @brief Gets the associated mesh to the simplex.
       */
      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

      // virtual VertexIterator getVertices() const;

      virtual SimplexIterator getAdjacent() const;

      virtual SimplexIterator getIncident() const;

      mfem::ElementTransformation& getTransformation() const;

      virtual std::vector<Geometry::Point> getIntegrationRule(int order) const;

    private:
      const size_t m_dimension;
      const Index m_index;
      std::reference_wrapper<const MeshBase> m_mesh;
      std::vector<Index> m_vertices;
      Attribute m_attr;
      Geometry::Type m_type;
      mutable std::unique_ptr<mfem::ElementTransformation> m_trans;
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
      Element(
          Index index,
          const MeshBase& mesh, const std::vector<Index>& vertices,
          Attribute attr = RODIN_DEFAULT_SIMPLEX_ATTRIBUTE);

      Element(const Element&) = delete;

      Element(Element&& other)
        :  Simplex(std::move(other))
      {}
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
      Face(
          Index index,
          const MeshBase& mesh, const std::vector<Index>& vertices,
          Attribute attr = RODIN_DEFAULT_SIMPLEX_ATTRIBUTE);

      Face(const Face&) = delete;

      Face(Face&& other)
        : Simplex(std::move(other))
      {}

      bool isBoundary() const;

      bool isInterface() const;
  };

  class Vertex : public Simplex
  {
    public:
      Vertex(
          Index index,
          const MeshBase& mesh,
          const Math::Vector& coordinates,
          Attribute attr = RODIN_DEFAULT_SIMPLEX_ATTRIBUTE);

      double x() const
      {
        assert(0 < m_coordinates.size());
        return operator()(0);
      }

      double y() const
      {
        assert(1 < m_coordinates.size());
        return operator()(1);
      }

      double z() const
      {
        assert(2 < m_coordinates.size());
        return operator()(2);
      }

      double operator()(size_t i) const;

      const Math::Vector& coordinates() const
      {
        return m_coordinates;
      }

    private:
      Math::Vector m_coordinates;
  };

  /**
   * @brief Represents a spatial point on a simplex.
   *
   * This class represents the point:
   * @f[
   *   p = x(r)
   * @f]
   * on some simplex @f$ \tau \in \mathcal{T}_h @f$ belonging to
   * some discrete mesh @f$ \mathcal{T}_h @f$. Here @f$ p \in \tau @f$ denotes
   * the physical coordinates of the point, while @f$ x : K \rightarrow \tau
   * @f$ represents the transformation taking reference coordinates @f$ r \in K
   * @f$, for a reference geometry @f$ K @f$.
   */
  class Point
  {
    public:
      enum class Coordinates
      {
        Reference,
        Physical
      };

      // Point(Coordinates type, const Transformation& trans, const Math::Vector& coords);

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
      size_t getDimension(Coordinates coords = Coordinates::Physical) const;

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

      // [[deprecated]]
      const mfem::IntegrationPoint& getIntegrationPoint() const
      {
        return m_ip;
      }

    private:
      mutable mfem::Vector m_physical;
      std::reference_wrapper<const Simplex> m_element;
      mfem::IntegrationPoint m_ip;
  };
}

#endif
