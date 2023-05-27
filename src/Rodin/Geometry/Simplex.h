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

#include "Rodin/Array.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"

#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class for all geometric elements of the mesh.
   */
  class Polytope
  {
    public:
      enum class Geometry
      {
        Point,
        Segment,
        Triangle,
        Quadrilateral,
        Tetrahedron
      };

      static constexpr std::array Geometries
      {
        Geometry::Point,
        Geometry::Segment,
        Geometry::Triangle,
        Geometry::Quadrilateral,
        Geometry::Tetrahedron
      };

      inline
      constexpr
      static size_t getVertexCount(Polytope::Geometry g)
      {
        switch (g)
        {
          case Geometry::Point:
            return 1;
          case Geometry::Segment:
            return 2;
          case Geometry::Triangle:
            return 3;
          case Geometry::Quadrilateral:
            return 4;
          case Geometry::Tetrahedron:
            return 4;
        }
        assert(false);
        return 0;
      }

      inline
      constexpr
      static size_t getGeometryDimension(Polytope::Geometry g)
      {
        switch (g)
        {
          case Geometry::Point:
            return 0;
          case Geometry::Segment:
            return 1;
          case Geometry::Triangle:
          case Geometry::Quadrilateral:
            return 2;
          case Geometry::Tetrahedron:
            return 3;
        }
        assert(false);
        return 0;
      }

      Polytope(size_t dimension, Index index, const MeshBase& mesh);

      Polytope(const Polytope&) = delete;

      Polytope(Polytope&&) = default;

      virtual ~Polytope() = default;

      /**
       * @brief Gets the index of the simplex in the mesh.
       */
      inline
      Index getIndex() const
      {
        return m_index;
      }

      inline
      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the associated mesh to the simplex.
       */
      inline
      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the attribute of the simplex.
       */
      Attribute getAttribute() const;

      Scalar getVolume() const;

      const PolytopeTransformation& getTransformation() const;

      VertexIterator getVertex() const;

      const Array<Index>& getVertices() const;

      PolytopeIterator getAdjacent() const;

      PolytopeIterator getIncident() const;

      Geometry getGeometry() const;

    private:
      const size_t m_dimension;
      const Index m_index;
      std::reference_wrapper<const MeshBase> m_mesh;
  };

  bool operator<(const Polytope& lhs, const Polytope& rhs);

  /**
   * @brief Class for representing elements of the highest dimension in the
   * mesh, i.e. tetrahedra in 3D, triangles in 2D or lines in 1D.
   *
   * This class is designed so that modifications cannot be made to the
   * element. If one wishes to modify the element then one must use
   * ElementView.
   */
  class Element : public Polytope
  {
    public:
      Element(Index index, const MeshBase& mesh);

      Element(const Element&) = delete;

      Element(Element&& other)
        :  Polytope(std::move(other))
      {}
  };

  /**
   * @brief Class for representing elements of codimension 1 in the
   * mesh, i.e. triangles in 3D or lines in 2D.
   *
   * This class is designed so that modifications cannot be made to the
   * face.
   */
  class Face : public Polytope
  {
    public:
      Face(Index index, const MeshBase& mesh);

      Face(const Face&) = delete;

      Face(Face&& other)
        : Polytope(std::move(other))
      {}

      bool isBoundary() const;

      bool isInterface() const;
  };

  class Vertex : public Polytope
  {
    public:
      Vertex(Index index, const MeshBase& mesh);

      inline
      Scalar x() const
      {
        return operator()(0);
      }

      inline
      Scalar y() const
      {
        return operator()(1);
      }

      inline
      Scalar z() const
      {
        return operator()(2);
      }

      inline
      Scalar operator()(size_t i) const
      {
        return getCoordinates()(i);
      }

      Eigen::Map<const Math::Vector> getCoordinates() const;

      inline
      constexpr
      Geometry getGeometry() const
      {
        return Geometry::Point;
      }
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
      Point(const Polytope& simplex, const PolytopeTransformation& trans, const Math::Vector& rc);

      Point(const Polytope& simplex, const PolytopeTransformation& trans, const Math::Vector& rc, const Math::Vector& pc);

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
            return getCoordinates(Coordinates::Reference)(i);
          }
        }
        return NAN;
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
      const Polytope& getSimplex() const
      {
        return m_simplex.get();
      }

      inline
      const PolytopeTransformation& getTransformation() const
      {
        return m_trans.get();
      }

      const Math::Vector& getCoordinates(Coordinates coords = Coordinates::Physical) const;

      const Math::Matrix& getJacobian() const;

      const Math::Matrix& getJacobianInverse() const;

      Scalar getDistortion() const;

    private:
      std::reference_wrapper<const Polytope> m_simplex;
      std::reference_wrapper<const PolytopeTransformation> m_trans;

      mutable std::optional<const Math::Vector> m_rc;
      mutable std::optional<const Math::Vector> m_pc;
      mutable std::optional<const Math::Matrix> m_jacobian;
      mutable std::optional<const Math::Matrix> m_inverseJacobian;
      mutable std::optional<const Scalar> m_distortion;
  };

  inline
  std::ostream& operator<<(std::ostream& os, Polytope::Geometry g)
  {
    switch (g)
    {
      case Polytope::Geometry::Point:
      {
        os << "Point";
        break;
      }
      case Polytope::Geometry::Segment:
      {
        os << "Segment";
        break;
      }
      case Polytope::Geometry::Triangle:
      {
        os << "Triangle";
        break;
      }
      case Polytope::Geometry::Quadrilateral:
      {
        os << "Quadrilateral";
        break;
      }
      case Polytope::Geometry::Tetrahedron:
      {
        os << "Tetrahedron";
        break;
      }
    }
    return os;
  }
}

#endif
