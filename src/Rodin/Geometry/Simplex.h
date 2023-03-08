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
#include <optional>
#include <mfem.hpp>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  enum class Type
  {
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
      inline
      Index getIndex() const
      {
        return m_index;
      }

      inline
      Type getGeometry() const
      {
        return m_type;
      }

      inline
      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the attribute of the simplex.
       */
      inline
      Attribute getAttribute() const
      {
        return m_attr;
      }

      /**
       * @brief Gets the associated mesh to the simplex.
       */
      inline
      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

      Scalar getVolume() const;

      const SimplexTransformation& getTransformation() const;

      // virtual VertexIterator getVertices() const;

      virtual SimplexIterator getAdjacent() const;

      virtual SimplexIterator getIncident() const;

    private:
      const size_t m_dimension;
      const Index m_index;
      std::reference_wrapper<const MeshBase> m_mesh;
      std::vector<Index> m_vertices;
      Attribute m_attr;
      Geometry::Type m_type;
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

      Scalar x() const
      {
        assert(0 < m_coordinates.size());
        return operator()(0);
      }

      Scalar y() const
      {
        assert(1 < m_coordinates.size());
        return operator()(1);
      }

      Scalar z() const
      {
        assert(2 < m_coordinates.size());
        return operator()(2);
      }

      Scalar operator()(size_t i) const;

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

      /**
       * @brief Constructs the Point object from reference coordinates.
       * @param[in] simplex Simplex to which point belongs to
       * @param[in] ip Reference coordinates
       */
      Point(const Simplex& simplex, const SimplexTransformation& trans, const Math::Vector& rc);

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
      inline
      Scalar operator()(size_t i, Coordinates coords = Coordinates::Physical) const
      {
        switch (coords)
        {
          case Coordinates::Physical:
          {
            return getCoordinates(Coordinates::Physical)(i);
          }
          case Coordinates::Reference:
          {
            assert(m_rc.get().size() > static_cast<int>(i));
            return m_rc.get()(i);
          }
        }
        return getCoordinates(Coordinates::Physical)(i);
      }

      /**
       * @brief Gets the @f$ x @f$ physical coordinate.
       * @returns Physical @f$ x @f$-coordinate.
       */
      inline
      Scalar x(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(0, coords);
      }

      /**
       * @brief Gets the @f$ y @f$ physical coordinate.
       * @returns Physical @f$ y @f$-coordinate.
       */
      inline
      Scalar y(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(1, coords);
      }

      /**
       * @brief Gets the @f$ z @f$ physical coordinate.
       * @returns Physical @f$ z @f$-coordinate.
       */
      inline
      Scalar z(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(2, coords);
      }

      /**
       * @brief Lexicographical comparison.
       */
      inline
      bool operator<(const Point& p) const
      {
        assert(getDimension() == p.getDimension());
        const Math::Vector& lhs = getCoordinates(Coordinates::Physical);
        const Math::Vector& rhs = p.getCoordinates(Coordinates::Physical);
        for (int i = 0; i < lhs.size() - 1; i++)
        {
          if (lhs(i) < rhs(i))
            return true;
          if (rhs(i) > lhs(i))
            return false;
        }
        return (lhs(lhs.size() - 1) < rhs(rhs.size() - 1));
      }

      inline
      const Simplex& getSimplex() const
      {
        return m_simplex.get();
      }

      inline
      const SimplexTransformation& getTransformation() const
      {
        return m_trans.get();
      }

      inline
      const mfem::IntegrationPoint& getIntegrationPoint() const
      {
        return m_ip;
      }

      const Math::Vector& getCoordinates(Coordinates coords = Coordinates::Physical) const;

      const Math::Matrix& getJacobian() const;

      const Math::Matrix& getJacobianInverse() const;

      Scalar getDistortion() const;

    private:
      std::reference_wrapper<const Simplex> m_simplex;
      std::reference_wrapper<const SimplexTransformation> m_trans;
      std::reference_wrapper<const Math::Vector> m_rc;
      mfem::IntegrationPoint m_ip;
      mutable std::optional<const Math::Vector> m_pc;
      mutable std::optional<const Math::Matrix> m_jacobian;
      mutable std::optional<const Math::Matrix> m_inverseJacobian;
      mutable std::optional<const Scalar> m_distortion;
  };
}

#endif
