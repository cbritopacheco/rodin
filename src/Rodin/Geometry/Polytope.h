/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPE_H
#define RODIN_GEOMETRY_POLYTOPE_H

#include <set>
#include <iostream>
#include <array>
#include <optional>

#include "Rodin/Configure.h"

#include "Rodin/Array.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Threads/Mutable.h"

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
        Tetrahedron,
        TriangularPrism
      };

      /**
       * @brief Iterable of possible polytope geometry types.
       */
      static constexpr std::array<Type, 6> Types
      {
        Type::Point,
        Type::Segment,
        Type::Triangle,
        Type::Quadrilateral,
        Type::Tetrahedron,
        Type::TriangularPrism
      };

      static const Math::PointMatrix& getVertices(Polytope::Type g);

      static auto getVertex(size_t i, Polytope::Type g)
      {
        return getVertices(g).col(i);
      }

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
          case Type::TriangularPrism:
            return 6;
        }
        assert(false);
        return 0;
      }

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
          case Type::TriangularPrism:
            return 3;
        }
        assert(false);
        return 0;
      }

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
          case Type::TriangularPrism:
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
      Index getIndex() const
      {
        return m_index;
      }

      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the associated mesh to the simplex.
       */
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
      Real getMeasure() const;

      const PolytopeTransformation& getTransformation() const;

      VertexIterator getVertex() const;

      const Array<Index>& getVertices() const;

      PolytopeIterator getAdjacent() const;

      Type getGeometry() const;

      bool isCell() const;

      bool isFace() const;

      bool isVertex() const;

      Polytope& setAttribute();

    private:
      static const GeometryIndexed<Math::PointMatrix> s_vertices;

      const size_t m_dimension;
      const Index m_index;
      std::reference_wrapper<const MeshBase> m_mesh;
  };

  bool operator==(const Polytope& lhs, const Polytope& rhs);

  bool operator<(const Polytope& lhs, const Polytope& rhs);

  /**
   * @brief Class for representing polytopes of the highest dimension in the
   * mesh, i.e. tetrahedra in 3D, triangles in 2D or lines in 1D.
   */
  class Cell : public Polytope
  {
    public:
      using Parent = Polytope;

      Cell(Index index, const MeshBase& mesh);

      Cell(const Cell& other)
        : Polytope(other)
      {}

      Cell(Cell&& other)
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

  /**
   * @brief Represents a vertex of the mesh.
   */
  class Vertex : public Polytope
  {
    public:
      using Parent = Polytope;

      Vertex(Index index, const MeshBase& mesh);

      Vertex(const Vertex& other)
        : Polytope(other)
      {}

      /**
       * @brief Acess the 1st-coordinate of the vertex.
       */
      Real x() const
      {
        return operator()(0);
      }

      /**
       * @brief Acess the 2nd-coordinate of the vertex.
       */
      Real y() const
      {
        return operator()(1);
      }

      /**
       * @brief Acess the 3rd-coordinate of the vertex.
       */
      Real z() const
      {
        return operator()(2);
      }

      /**
       * @brief Acess the ith-coordinate of the vertex.
       */
      Real operator()(size_t i) const
      {
        return getCoordinates()(i);
      }

      Eigen::Map<const Math::SpatialVector<Real>> getCoordinates() const;

      constexpr
      Type getGeometry() const
      {
        return Type::Point;
      }
  };

  std::ostream& operator<<(std::ostream& os, const Polytope::Type& p);
}

#endif
