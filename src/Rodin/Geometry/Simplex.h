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
       * @brief The type of the Polytope geometry.
       */
      enum class Type
      {
        Point,
        Segment,
        Triangle,
        Quadrilateral,
        Tetrahedron
      };

      /**
       * @brief Iterable of possible polytope geometry types.
       */
      static constexpr std::array Types
      {
        Type::Point,
        Type::Segment,
        Type::Triangle,
        Type::Quadrilateral,
        Type::Tetrahedron
      };

      static const Math::Matrix& getVertices(Polytope::Type g);

      inline
      constexpr
      static size_t getVertexCount(Polytope::Type g)
      {
        switch (g)
        {
          case Type::Point:
            return 1;
          case Type::Segment:
            return 2;
          case Type::Triangle:
            return 3;
          case Type::Quadrilateral:
          case Type::Tetrahedron:
            return 4;
        }
        assert(false);
        return 0;
      }

      inline
      constexpr
      static size_t getGeometryDimension(Polytope::Type g)
      {
        switch (g)
        {
          case Type::Point:
            return 0;
          case Type::Segment:
            return 1;
          case Type::Triangle:
          case Type::Quadrilateral:
            return 2;
          case Type::Tetrahedron:
            return 3;
        }
        assert(false);
        return 0;
      }

      inline
      constexpr
      static bool isSimplex(Polytope::Type g)
      {
        switch (g)
        {
          case Type::Point:
          case Type::Segment:
          case Type::Triangle:
          case Type::Tetrahedron:
            return true;
          case Type::Quadrilateral:
            return false;
        }
        assert(false);
        return false;
      }

      /**
       * @brief Consructs a polytope of dimension @f$ d @f$ and index @f$ i @f$
       * belonging to the given mesh.
       */
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

      /**
       * @brief Gets the measure of the polytope.
       *
       * Gets the @f$ d @f$-dimensional measure of the polytope. This has
       * different names in different dimensions. See table below.
       * Dimension of polytope    | Measure
       * ------------------------ | -------------
       * 0                        | Always zero
       * 1                        | Length
       * 2                        | Area
       * 3                        | Volume
       *
       * @return The measure of the polytope.
       */
      Scalar getMeasure() const;

      const PolytopeTransformation& getTransformation() const;

      VertexIterator getVertex() const;

      const Array<Index>& getVertices() const;

      PolytopeIterator getAdjacent() const;

      PolytopeIterator getIncident() const;

      Type getGeometry() const;

    private:
      static const GeometryIndexed<Math::Matrix> s_vertices;

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

      Eigen::Map<const Math::SpatialVector> getCoordinates() const;

      inline
      constexpr
      Type getGeometry() const
      {
        return Type::Point;
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
    enum class PolytopeStorage
    {
      Value,
      Reference
    };

    public:
      /// Denotes the type of coordinates.
      enum class Coordinates
      {
        Reference, ///< Reference coordinates
        Physical ///< Physical coordinates
      };

      explicit
      PointBase(std::reference_wrapper<const Polytope> polytope, const PolytopeTransformation& trans,
          const Math::SpatialVector& pc);

      explicit
      PointBase(std::reference_wrapper<const Polytope> polytope, const PolytopeTransformation& trans);

      explicit
      PointBase(Polytope&& polytope, const PolytopeTransformation& trans,
          const Math::SpatialVector& pc);

      explicit
      PointBase(Polytope&& polytope, const PolytopeTransformation& trans);

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
      bool operator<(const PointBase& p) const;

      const Polytope& getPolytope() const;

      inline
      const PolytopeTransformation& getTransformation() const
      {
        return m_trans.get();
      }

      const Math::SpatialVector& getPhysicalCoordinates() const;

      const Math::SpatialVector& getCoordinates(Coordinates coords = Coordinates::Physical) const;

      /**
       * @brief Computes the Jacobian matrix of the transformation at the
       * point.
       */
      virtual const Math::SpatialMatrix& getJacobian() const;

      Scalar getJacobianDeterminant() const;

      /**
       * @brief Computes the inverse of the Jacobian matrix of the
       * transformation at the point.
       */
      const Math::SpatialMatrix& getJacobianInverse() const;

      /**
       * @brief Computes the distortion of space of the transformation at the
       * point.
       */
      Scalar getDistortion() const;

      virtual const Math::SpatialVector& getReferenceCoordinates() const = 0;

    private:
      PolytopeStorage m_polytopeStorage;
      std::variant<const Polytope, std::reference_wrapper<const Polytope>> m_polytope;
      std::reference_wrapper<const PolytopeTransformation> m_trans;

      mutable std::optional<const Math::SpatialVector> m_pc;
      mutable std::optional<const Math::SpatialMatrix> m_jacobian;
      mutable std::optional<const Math::SpatialMatrix> m_jacobianInverse;
      mutable std::optional<const Scalar>              m_jacobianDeterminant;
      mutable std::optional<const Scalar>              m_distortion;
  };

  class Point final : public PointBase
  {
    enum class RCStorage
    {
      Value,
      Reference
    };

    public:
      using Parent = PointBase;

      explicit
      Point(std::reference_wrapper<const Polytope> polytope, const PolytopeTransformation& trans,
          std::reference_wrapper<const Math::SpatialVector> rc, const Math::SpatialVector& pc);

      explicit
      Point(std::reference_wrapper<const Polytope> polytope, const PolytopeTransformation& trans,
          Math::SpatialVector&& rc, const Math::SpatialVector& pc);

      explicit
      Point(std::reference_wrapper<const Polytope> polytope, const PolytopeTransformation& trans,
          std::reference_wrapper<const Math::SpatialVector> rc);

      explicit
      Point(std::reference_wrapper<const Polytope> polytope, const PolytopeTransformation& trans,
          Math::SpatialVector&& rc);

      explicit
      Point(Polytope&& polytope, const PolytopeTransformation& trans,
          std::reference_wrapper<const Math::SpatialVector> rc, const Math::SpatialVector& pc);

      explicit
      Point(Polytope&& polytope, const PolytopeTransformation& trans,
          Math::SpatialVector&& rc, const Math::SpatialVector& pc);

      explicit
      Point(Polytope&& polytope, const PolytopeTransformation& trans,
          std::reference_wrapper<const Math::SpatialVector> rc);

      explicit
      Point(Polytope&& polytope, const PolytopeTransformation& trans,
          Math::SpatialVector&& rc);

      Point(const Point& other);

      Point(Point&& other);

      const Math::SpatialVector& getReferenceCoordinates() const override;

    private:
      const RCStorage m_rcStorage;
      std::variant<const Math::SpatialVector, std::reference_wrapper<const Math::SpatialVector>> m_rc;
  };

  inline
  std::ostream& operator<<(std::ostream& os, Polytope::Type g)
  {
    switch (g)
    {
      case Polytope::Type::Point:
      {
        os << "Point";
        break;
      }
      case Polytope::Type::Segment:
      {
        os << "Segment";
        break;
      }
      case Polytope::Type::Triangle:
      {
        os << "Triangle";
        break;
      }
      case Polytope::Type::Quadrilateral:
      {
        os << "Quadrilateral";
        break;
      }
      case Polytope::Type::Tetrahedron:
      {
        os << "Tetrahedron";
        break;
      }
    }
    return os;
  }
}

#endif
