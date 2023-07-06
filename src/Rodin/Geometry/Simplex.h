/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEX_H
#define RODIN_GEOMETRY_SIMPLEX_H

#include <set>
#include <iostream>
#include <array>
#include <optional>

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
      /**
       * @brief Polytope geometry
       */
      enum class Geometry
      {
        Point,
        Segment,
        Triangle,
        Quadrilateral,
        Tetrahedron
      };

      /**
       * @brief Iterable of possible polytope geometries.
       */
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

      inline
      constexpr
      static bool isSimplex(Polytope::Geometry g)
      {
        switch (g)
        {
          case Geometry::Point:
          case Geometry::Segment:
          case Geometry::Triangle:
          case Geometry::Tetrahedron:
            return true;
          default:
            return false;
        }
        assert(false);
        return false;
      }

      Polytope(size_t dimension, Index index, const MeshBase& mesh);

      Polytope(const Polytope&) = default;

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
      using Parent = Polytope;

      Element(Index index, const MeshBase& mesh);

      Element(const Element& other)
        : Polytope(other)
      {}

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
      using Parent = Polytope;

      Face(Index index, const MeshBase& mesh);

      Face(const Face& other)
        : Polytope(other)
      {}

      Face(Face&& other)
        : Polytope(std::move(other))
      {}

      bool isBoundary() const;

      bool isInterface() const;
  };

  class Vertex : public Polytope
  {
    public:
      using Parent = Polytope;

      Vertex(Index index, const MeshBase& mesh);

      Vertex(const Vertex& other)
        : Polytope(other)
      {}

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
   * @brief Abstract class for spatial points on a discrete mesh.
   *
   * This class represents the tuple @f$ (x, r, p) @f$
   * such that:
   * @f[
   *  p = x(r)
   * @f]
   * for a polytope @f$ \tau \in \mathcal{T}_h @f$ belonging to the mesh @f$
   * \mathcal{T}_h @f$. Here, @f$ p \in \tau @f$ denotes the physical
   * coordinates of the point, while @f$ x : K \rightarrow \tau @f$ represents
   * the transformation taking reference coordinates @f$ r \in K @f$, for a
   * reference geometry @f$ K @f$.
   *
   * @section rodin-geometry-point-thread_safety Thread safety
   * This class is not thread safe.
   *
   * @see PolytopeTransformation
   */
  class PointBase
  {
    public:
      /// Denotes the type of coordinates.
      enum class Coordinates
      {
        Reference, ///< Reference coordinates
        Physical ///< Physical coordinates
      };

      PointBase(const Polytope& polytope, const PolytopeTransformation& trans);

      PointBase(const Polytope& polytope, const PolytopeTransformation& trans, const Math::Vector& pc);

      /**
       * @brief Copy constructor.
       */
      PointBase(const PointBase&) = default;

      /**
       * @brief Move constructor.
       */
      PointBase(PointBase&&) = default;

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
            return getPhysicalCoordinates()(i);
          }
          case Coordinates::Reference:
          {
            return getReferenceCoordinates()(i);
          }
        }
        return NAN;
      }

      /**
       * @brief Gets the @f$ x @f$ coordinate.
       * @returns @f$ x @f$ coordinate of the point.
       */
      inline
      Scalar x(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(0, coords);
      }

      /**
       * @brief Gets the @f$ y @f$ coordinate.
       * @returns @f$ y @f$ coordinate of the point.
       */
      inline
      Scalar y(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(1, coords);
      }

      /**
       * @brief Gets the @f$ z @f$ coordinate.
       * @returns @f$ z @f$ coordinate of the point.
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
      bool operator<(const PointBase& p) const
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
      const Polytope& getPolytope() const
      {
        return m_polytope.get();
      }

      inline
      const PolytopeTransformation& getTransformation() const
      {
        return m_trans.get();
      }

      const Math::Vector& getPhysicalCoordinates() const;

      const Math::Vector& getCoordinates(Coordinates coords = Coordinates::Physical) const;

      /**
       * @brief Computes the Jacobian matrix of the transformation at the
       * point.
       */
      virtual const Math::Matrix& getJacobian() const;

      Scalar getJacobianDeterminant() const;

      /**
       * @brief Computes the inverse of the Jacobian matrix of the
       * transformation at the point.
       */
      const Math::Matrix& getJacobianInverse() const;

      /**
       * @brief Computes the distortion of space of the transformation at the
       * point.
       */
      Scalar getDistortion() const;

      virtual const Math::Vector& getReferenceCoordinates() const = 0;

    private:
      std::reference_wrapper<const Polytope> m_polytope;
      std::reference_wrapper<const PolytopeTransformation> m_trans;

      mutable std::optional<const Math::Vector> m_pc;
      mutable std::optional<const Math::Matrix> m_jacobian;
      mutable std::optional<const Math::Matrix> m_jacobianInverse;
      mutable std::optional<const Scalar>       m_jacobianDeterminant;
      mutable std::optional<const Scalar>       m_distortion;
  };

  class Point final : public PointBase
  {
    public:
      using Parent = PointBase;

      enum class Type
      {
        Data,
        Reference
      };

      Point(const Polytope& polytope, const PolytopeTransformation& trans, const Math::Vector& rc)
        : PointBase(polytope, trans), m_type(Type::Data), m_rc(rc)
      {}

      explicit
      Point(const Polytope& polytope, const PolytopeTransformation& trans, std::reference_wrapper<const Math::Vector> rc)
        : PointBase(polytope, trans), m_type(Type::Reference), m_rc(rc)
      {}

      Point(const Point& other)
        : PointBase(other),
          m_type(other.m_type),
          m_rc(other.m_rc)
      {}

      Point(Point&& other)
        : PointBase(std::move(other)),
          m_type(other.m_type),
          m_rc(std::move(other.m_rc))
      {}

      inline
      constexpr
      bool holds(Type t) const
      {
        return m_type == t;
      }

      inline
      const Math::Vector& getReferenceCoordinates() const override
      {
        if (holds(Type::Data))
        {
          return std::get<const Math::Vector>(m_rc);
        }
        else
        {
          assert(holds(Type::Reference));
          return std::get<std::reference_wrapper<const Math::Vector>>(m_rc);
        }
      }

    private:
      const Type m_type;
      std::variant<const Math::Vector, std::reference_wrapper<const Math::Vector>> m_rc;
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
