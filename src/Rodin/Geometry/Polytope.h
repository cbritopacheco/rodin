/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPE_H
#define RODIN_GEOMETRY_POLYTOPE_H

/**
 * @file
 * @brief Polytope classes representing geometric elements in finite element meshes.
 */

#include <iostream>
#include <array>

#include "Rodin/Configure.h"

#include "Rodin/Array.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

#include "ForwardDecls.h"

#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class for all geometric elements in finite element meshes.
   *
   * A Polytope represents a geometric entity in a finite element mesh, ranging
   * from vertices (0D) to cells (highest dimension). Polytopes are the building
   * blocks of mesh topology and provide the geometric foundation for finite
   * element computations.
   *
   * # Mathematical Foundation
   *
   * Each polytope @f$ \tau @f$ in a mesh @f$ \mathcal{T}_h @f$ has:
   * - A dimension @f$ d @f$ (0 for vertices, 1 for edges, etc.)
   * - An index @f$ i @f$ unique within its dimension
   * - An associated reference element @f$ K @f$
   * - A transformation @f$ x: K \rightarrow \tau @f$ mapping reference to physical coordinates
   *
   * # Hierarchy
   *
   * The Polytope class serves as the base for specialized classes:
   * - Cell: Highest-dimensional elements (e.g., tetrahedra in 3D, triangles in 2D)
   * - Face: Codimension-1 elements (e.g., triangles in 3D, edges in 2D)
   * - Vertex: 0-dimensional points
   *
   * # Usage
   *
   * Polytopes are typically obtained through iterators:
   * @code{.cpp}
   * for (auto it = mesh.getCell(); it; ++it)
   * {
   *   const Cell& cell = *it;
   *   std::cout << "Cell " << cell.getIndex() 
   *             << " has measure " << cell.getMeasure() << std::endl;
   * }
   * @endcode
   *S 
   * @see Cell, Face, Vertex, PolytopeIterator
   */
  class Polytope
  {
    public:
      class Key
      {
        public:
          using Vertices = std::array<Index, RODIN_MAXIMUM_POLYTOPE_VERTICES>;
          static_assert(RODIN_MAXIMUM_POLYTOPE_VERTICES == 8);

          using iterator = Vertices::iterator;
          using const_iterator = Vertices::const_iterator;
          using value_type = Index;
          using reference = Index&;
          using const_reference = const Index&;
          using size_type = std::uint8_t;

          struct SymmetricEquality
          {
            bool operator()(const Key& a, const Key& b) const;
          };

          struct SymmetricHash
          {
            std::size_t operator()(const Key& poly) const;
          };

          constexpr
          static inline std::uint64_t sm64(std::uint64_t x)
          {
            x += 0x9e3779b97f4a7c15ull;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
            x ^= (x >> 31);
            return x;
          }

          constexpr
          static inline void cswap(Index& a, Index& b)
          {
            if (b < a)
              std::swap(a, b);
          }

          Key();

          Key(std::uint8_t n);

          Key(std::initializer_list<Index> vertices);

          Index& operator()(std::uint8_t i);

          const Index& operator()(std::uint8_t i) const;

          Index& operator[](std::uint8_t i);

          const Index& operator[](std::uint8_t i) const;

          Key& resize(std::uint8_t n);

          std::uint8_t size() const;

          const Vertices& getVertices() const;

          Vertices::iterator begin();

          Vertices::iterator end();

          Vertices::const_iterator begin() const;

          Vertices::const_iterator end() const;

          template <class Archive>
          void serialize(Archive& ar, const unsigned int)
          {
            ar & m_n;
            ar & m_vertices;
          }

        private:
          std::uint8_t m_n;
          Vertices m_vertices;
      };

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
        Hexahedron,     ///< 3D hexahedral element
        Wedge           ///< 3D prismatic element (triangular prism)
      };

      /**
       * @brief Projection operations for polytope sub-entities.
       *
       * Provides methods to project points from one polytope's reference
       * coordinates to another (e.g., from a face to its parent cell).
       */
      struct Project
      {
        public:
          /**
           * @brief Constructs a projector for a given geometry type.
           * @param[in] g Geometry type
           */
          Project(Type g);

          /**
           * @brief Projects reference coordinates to cell coordinates.
           * @param[out] out Output coordinates
           * @param[in] rc Reference coordinates
           */
          void cell(Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;

          /**
           * @brief Projects reference coordinates to boundary coordinates.
           * @param[out] out Output coordinates
           * @param[in] rc Reference coordinates
           */
          void boundary(Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;

          /**
           * @brief Projects to a specific face.
           * @param[in] local Local face index
           * @param[out] out Output coordinates
           * @param[in] rc Reference coordinates
           */
          void face(size_t local, Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;

          /**
           * @brief Projects to a specific vertex.
           * @param[in] local Local vertex index
           * @param[out] out Output coordinates
           * @param[in] rc Reference coordinates
           */
          void vertex(size_t local, Math::SpatialPoint& out, const Math::SpatialPoint& rc) const;
        private:
          const Type m_g;
      };

      /**
       * @brief Geometric traits for polytope types.
       *
       * Provides compile-time and runtime information about polytope geometries,
       * including dimension, vertex count, and reference element properties.
       */
      struct Traits
      {
        public:
          /**
           * @brief Half-space representation of reference element.
           *
           * Defines the reference element using linear inequalities
           * @f$ Ax \leq b @f$.
           */
          struct HalfSpace
          {
            Math::Matrix<Real> matrix; ///< Matrix A in Ax ≤ b
            Math::Vector<Real> vector; ///< Vector b in Ax ≤ b
          };

          /**
           * @brief Constructs traits for a geometry type.
           * @param[in] g Geometry type
           */
          Traits(Type g);

          /**
           * @brief Checks if the geometry is a simplex.
           * @returns True for triangles and tetrahedra, false otherwise
           */
          bool isSimplex() const;

          bool isTensorProduct() const;

          /**
           * @brief Gets the topological dimension.
           * @returns Dimension (0 for points, 1 for segments, 2 for triangles/quads, 3 for tets/wedges)
           */
          size_t getDimension() const;

          /**
           * @brief Gets the number of vertices.
           * @returns Vertex count for this geometry type
           */
          size_t getVertexCount() const;

          /**
           * @brief Gets reference coordinates of a vertex.
           * @param[in] i Vertex index
           * @returns Reference coordinates of vertex i
           */
          const Math::SpatialPoint& getVertex(size_t i) const;

          /**
           * @brief Gets the half-space representation.
           * @returns Half-space defining the reference element
           */
          const HalfSpace& getHalfSpace() const;

          Math::SpatialPoint getCentroid() const;

        private:
          const Type m_g;
      };

      /**
       * @brief Array of all supported polytope types.
       *
       * Useful for iterating over all geometry types at compile time.
       */
      static constexpr std::array<Type, 7> Types
      {
        Type::Point,
        Type::Segment,
        Type::Triangle,
        Type::Quadrilateral,
        Type::Tetrahedron,
        Type::Hexahedron,
        Type::Wedge
      };

      /**
       * @brief Constructs a polytope with given dimension and index.
       * @param[in] dimension Topological dimension @f$ d @f$
       * @param[in] index Index @f$ i @f$ within the mesh
       * @param[in] mesh Reference to the containing mesh
       */
      Polytope(size_t dimension, Index index, const MeshBase& mesh)
        : m_dimension(dimension), m_index(index), m_mesh(mesh)
      {}

      /**
       * @brief Copy constructor.
       */
      Polytope(const Polytope& other)
        : m_dimension(other.m_dimension),
          m_index(other.m_index),
          m_mesh(other.m_mesh)
      {}

      /**
       * @brief Move constructor.
       */
      Polytope(Polytope&& other)
        : m_dimension(std::exchange(other.m_dimension, 0)),
          m_index(std::exchange(other.m_index, 0)),
          m_mesh(std::move(other.m_mesh))
      {}

      /**
       * @brief Copy assignment operator.
       */
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

      /**
       * @brief Move assignment operator.
       */
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

      /**
       * @brief Virtual destructor.
       */
      virtual ~Polytope() = default;

      /**
       * @brief Gets the polytope's index within its dimension.
       * @returns Index @f$ i @f$ of this polytope
       */
      Index getIndex() const
      {
        return m_index;
      }

      /**
       * @brief Gets the topological dimension.
       * @returns Dimension @f$ d @f$ of this polytope
       */
      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the containing mesh.
       * @returns Reference to the mesh this polytope belongs to
       */
      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the attribute (marker) of this polytope.
       * @returns Attribute value
       *
       * Attributes are typically used to mark material regions, boundary
       * conditions, or other domain-specific properties.
       */
      Optional<Attribute> getAttribute() const;

      /**
       * @brief Gets the geometric measure of this polytope.
       *
       * Computes the @f$ d @f$-dimensional measure of the polytope:
       *
       * | Dimension | Measure    |
       * |-----------|------------|
       * | 0         | Always 0   |
       * | 1         | Length     |
       * | 2         | Area       |
       * | 3         | Volume     |
       *
       * @returns The measure (length/area/volume) of the polytope
       */
      Real getMeasure() const;

      /**
       * @brief Gets the geometric transformation for this polytope.
       * @returns Reference to the transformation @f$ x: K \rightarrow \tau @f$
       *
       * @see PolytopeTransformation
       */
      const PolytopeTransformation& getTransformation() const;

      /**
       * @brief Gets an iterator over the vertices of this polytope.
       * @returns Iterator over all vertices
       */
      VertexIterator getVertex() const;

      /**
       * @brief Gets the vertex indices defining this polytope.
       * @returns Array of vertex indices
       */
      const Key& getVertices() const;

      /**
       * @brief Gets an iterator over polytopes adjacent to this one.
       * @returns Iterator over adjacent polytopes
       *
       * @note Requires connectivity to be computed first.
       */
      PolytopeIterator getAdjacent() const;

      /**
       * @brief Gets the geometry type of this polytope.
       * @returns Geometry type (Triangle, Tetrahedron, etc.)
       */
      Type getGeometry() const;

      /**
       * @brief Checks if this is a cell (maximal dimension).
       * @returns True if this polytope is a cell, false otherwise
       */
      bool isCell() const;

      /**
       * @brief Checks if this is a face (codimension 1).
       * @returns True if this polytope is a face, false otherwise
       */
      bool isFace() const;

      /**
       * @brief Checks if this is a vertex (dimension 0).
       * @returns True if this polytope is a vertex, false otherwise
       */
      bool isVertex() const;

      /**
       * @brief Sets the attribute for this polytope.
       * @returns Reference to this polytope for method chaining
       *
       * @note Implementation details in derived classes.
       */
      Polytope& setAttribute();

    private:
      size_t m_dimension;
      Index m_index;
      std::reference_wrapper<const MeshBase> m_mesh;
  };

  /**
   * @brief Equality comparison for polytopes.
   * @param[in] lhs Left-hand side polytope
   * @param[in] rhs Right-hand side polytope
   * @returns True if polytopes are equal (same mesh, dimension, and index)
   */
  bool operator==(const Polytope& lhs, const Polytope& rhs);

  /**
   * @brief Less-than comparison for polytopes.
   * @param[in] lhs Left-hand side polytope
   * @param[in] rhs Right-hand side polytope
   * @returns True if lhs < rhs (lexicographic order by dimension then index)
   */
  bool operator<(const Polytope& lhs, const Polytope& rhs);

  /**
   * @brief Represents a cell (highest-dimensional element) in a mesh.
   *
   * Cells are polytopes of maximal dimension in the mesh:
   * - Tetrahedra in 3D meshes
   * - Triangles or quadrilaterals in 2D meshes  
   * - Segments in 1D meshes
   *
   * Cells form the basis for finite element computations and define the
   * computational domain.
   *
   * @see Polytope, Face, Vertex
   */
  class Cell : public Polytope
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = Polytope;

      /**
       * @brief Constructs a cell with given index.
       * @param[in] index Index of the cell in the mesh
       * @param[in] mesh Reference to the containing mesh
       */
      Cell(Index index, const MeshBase& mesh);

      /**
       * @brief Copy constructor.
       */
      Cell(const Cell& other)
        : Polytope(other)
      {}

      /**
       * @brief Move constructor.
       */
      Cell(Cell&& other)
        :  Polytope(std::move(other))
      {}
  };

  /**
   * @brief Represents a face (codimension-1 element) in a mesh.
   *
   * Faces are polytopes of dimension @f$ d-1 @f$ where @f$ d @f$ is the
   * mesh dimension:
   * - Triangles or quadrilaterals in 3D meshes
   * - Segments in 2D meshes
   *
   * Faces are crucial for:
   * - Boundary condition application
   * - Interface problems
   * - Flux computations
   *
   * @see Polytope, Cell, Vertex
   */
  class Face : public Polytope
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = Polytope;

      /**
       * @brief Constructs a face with given index.
       * @param[in] index Index of the face in the mesh
       * @param[in] mesh Reference to the containing mesh
       */
      Face(Index index, const MeshBase& mesh);

      /**
       * @brief Copy constructor.
       */
      Face(const Face& other)
        : Polytope(other)
      {}

      /**
       * @brief Move constructor.
       */
      Face(Face&& other)
        : Polytope(std::move(other))
      {}

      /**
       * @brief Checks if this face is on the mesh boundary.
       * @returns True if the face is a boundary face, false otherwise
       *
       * A face is on the boundary if it belongs to only one cell.
       */
      bool isBoundary() const;

      /**
       * @brief Checks if this face is an interface between regions.
       * @returns True if the face separates different material regions
       *
       * An interface face connects cells with different attributes.
       */
      bool isInterface() const;
  };

  /**
   * @brief Represents a vertex (0-dimensional element) in a mesh.
   *
   * Vertices are the fundamental building blocks of the mesh, representing
   * points in space with associated coordinates. They serve as:
   * - Degrees of freedom locations for P1 finite elements
   * - Endpoints for edges and corners for higher-dimensional elements
   * - Spatial anchors for the mesh geometry
   *
   * @see Polytope, Cell, Face
   */
  class Vertex : public Polytope
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = Polytope;

      /**
       * @brief Constructs a vertex with given index.
       * @param[in] index Index of the vertex in the mesh
       * @param[in] mesh Reference to the containing mesh
       */
      Vertex(Index index, const MeshBase& mesh);

      /**
       * @brief Copy constructor.
       */
      Vertex(const Vertex& other)
        : Polytope(other)
      {}

      /**
       * @brief Gets the x-coordinate (first coordinate).
       * @returns x-coordinate value
       */
      Real x() const
      {
        return operator()(0);
      }

      /**
       * @brief Gets the y-coordinate (second coordinate).
       * @returns y-coordinate value
       */
      Real y() const
      {
        return operator()(1);
      }

      /**
       * @brief Gets the z-coordinate (third coordinate).
       * @returns z-coordinate value
       */
      Real z() const
      {
        return operator()(2);
      }

      /**
       * @brief Gets the i-th coordinate.
       * @param[in] i Coordinate index
       * @returns i-th coordinate value
       */
      Real operator()(size_t i) const
      {
        return getCoordinates()(i);
      }

      /**
       * @brief Gets the vertex coordinates as a vector.
       * @returns Eigen map to the coordinate vector
       */
      Math::SpatialPoint getCoordinates() const;

      /**
       * @brief Gets the geometry type (always Point for vertices).
       * @returns Polytope::Type::Point
       */
      constexpr
      Type getGeometry() const
      {
        return Type::Point;
      }
  };

  /**
   * @brief Stream output operator for polytope types.
   * @param[in,out] os Output stream
   * @param[in] p Polytope type to output
   * @returns Reference to the output stream
   */
  std::ostream& operator<<(std::ostream& os, const Polytope::Type& p);
}

#endif
