/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPE_H
#define RODIN_GEOMETRY_POLYTOPE_H

#include <iostream>
#include <array>

#include "Rodin/Configure.h"

#include "Rodin/Array.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"

#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @ingroup RodinGeometry
   * @brief Base class for all geometric elements in finite element meshes.
   *
   * A Polytope represents a geometric entity in a finite element mesh, ranging
   * from vertices (0D) to cells (highest dimension). Polytopes are the building
   * blocks of mesh topology and provide the geometric foundation for finite
   * element computations.
   */
  class Polytope
  {
    public:
      /**
       * @brief Enumeration of supported polytope geometries.
       *
       * This enumeration defines the geometric types supported in Rodin's
       * finite element framework, covering both simplicial and tensor-product
       * element families commonly used in numerical analysis.
       */
      enum class Type
      {
        Point,          ///< 0D vertex element
        Segment,        ///< 1D line element
        Triangle,       ///< 2D triangular element
        Quadrilateral,  ///< 2D quadrilateral element
        Tetrahedron,    ///< 3D tetrahedral element
        Wedge           ///< 3D prismatic element
      };

      struct Project
      {
        public:
          Project(Type g);
          void cell(Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;
          void boundary(Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;
          void face(size_t local, Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;
          void vertex(size_t local, Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;
        private:
          const Type m_g;
      };

      struct Traits
      {
        public:
          struct HalfSpace
          {
            Math::Matrix<Real> matrix;
            Math::Vector<Real> vector;
          };

          Traits(Type g);

          bool isSimplex() const;

          size_t getDimension() const;

          size_t getVertexCount() const;

          const Math::SpatialPoint& getVertex(size_t i) const;

          const HalfSpace& getHalfSpace() const;

        private:
          const Type m_g;
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
        Type::Wedge
      };

      /**
       * @brief Consructs a polytope of dimension @f$ d @f$ and index @f$ i @f$
       * belonging to the given mesh.
       */
      Polytope(size_t dimension, Index index, const MeshBase& mesh)
        : m_dimension(dimension), m_index(index), m_mesh(mesh)
      {}

      Polytope(const Polytope& other)
        : m_dimension(other.m_dimension),
          m_index(other.m_index),
          m_mesh(other.m_mesh)
      {}

      Polytope(Polytope&& other)
        : m_dimension(std::exchange(other.m_dimension, 0)),
          m_index(std::exchange(other.m_index, 0)),
          m_mesh(std::move(other.m_mesh))
      {}

      Polytope& operator=(const Polytope& other)
      {
        if (this != &other)
        {
          m_dimension = other.m_dimension;
          m_index = other.m_index;
          m_mesh = other.m_mesh;
        }
        return *this;
      }

      Polytope& operator=(Polytope&& other)
      {
        if (this != &other)
        {
          m_dimension = std::exchange(other.m_dimension, 0);
          m_index = std::exchange(other.m_index, 0);
          m_mesh = std::move(other.m_mesh);
        }
        return *this;
      }

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
      size_t m_dimension;
      Index m_index;
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

      Eigen::Map<const Math::SpatialPoint> getCoordinates() const;

      constexpr
      Type getGeometry() const
      {
        return Type::Point;
      }
  };

  std::ostream& operator<<(std::ostream& os, const Polytope::Type& p);
}

#endif
