/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MESH_H
#define RODIN_GEOMETRY_MESH_H

/**
 * @file
 * @brief Core mesh data structure for finite element computations.
 *
 * This file defines the main Mesh class and related utilities for representing
 * and manipulating unstructured finite element meshes. The Mesh class provides:
 * - Topology management through connectivity relations
 * - Geometry storage and access
 * - Mesh generation and I/O operations
 * - Iterators for traversing mesh entities
 * - Support for submeshes and domain decomposition
 * - Attribute-based region marking
 */

#include <deque>

#include <boost/filesystem.hpp>
#include <boost/serialization/access.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Context/Local.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Geometry/AttributeIndex.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Rodin/Serialization/Array.h"
#include "Rodin/Serialization/EigenMatrix.h"

#include "Rodin/Alert/MemberFunctionException.h"

#include "ForwardDecls.h"
#include "Connectivity.h"
#include "Point.h"
#include "Polytope.h"
#include "PolytopeIterator.h"
#include "PolytopeTransformation.h"
#include "PolytopeTransformationIndex.h"

/**
 * @ingroup RodinDirectives
 * @brief Requires connectivity @f$ d \rightarrow d' @f$ to be computed.
 *
 * Throws an exception if the specified connectivity relation has not been
 * computed for the mesh. Connectivity must be explicitly computed using
 * `mesh.getConnectivity().compute(d, dp)` before operations that depend on it.
 *
 * @param[in] mesh Mesh instance to check
 * @param[in] d Source dimension
 * @param[in] dp Target dimension
 */
#define RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, dp) \
  if (mesh.getConnectivity().getIncidence(d, dp).size() == 0) \
  { \
    Rodin::Alert::MemberFunctionException(*this, __func__) \
      << Rodin::Alert::Notation::Incidence(d, dp) \
      << " has not been computed and is required to use this function." \
      << Rodin::Alert::Raise; \
  }

/**
 * @ingroup RodinDirectives
 * @brief Requires connectivity @f$ d \rightarrow d' @f$ (for use in member functions).
 *
 * Throws an exception if the specified connectivity relation has not been
 * computed. This macro is for use inside mesh member functions.
 *
 * @param[in] d Source dimension
 * @param[in] dp Target dimension
 */
#define RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(d, dp) \
  if (this->getConnectivity().getIncidence(d, dp).size() == 0) \
  { \
    Rodin::Alert::MemberFunctionException(*this, __func__) \
      << Rodin::Alert::Notation::Incidence(d, dp) \
      << " has not been computed and is required to use this function." \
      << Rodin::Alert::Raise; \
  }

/**
 * @ingroup RodinDirectives
 * @brief Requires mesh to be a SubMesh.
 *
 * Throws an exception if the mesh is not a SubMesh instance, preventing
 * invalid downcasting operations.
 */
#define RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH() \
  if (!this->isSubMesh()) \
  { \
    Rodin::Alert::MemberFunctionException(*this, __func__) \
      << "This instance of Mesh is not a SubMesh " \
      << Rodin::Alert::Notation::Predicate(false, "isSubMesh()") \
      << ". Downcasting to SubMesh is ill-defined." \
      << Rodin::Alert::Raise; \
  }

namespace Rodin::FormLanguage
{
  /**
   * @brief Traits specialization for Mesh types.
   * @tparam Context Execution context (e.g., Context::Local, Context::MPI)
   *
   * Provides compile-time information about mesh type properties for the
   * form language system.
   */
  template <class Context>
  struct Traits<Geometry::Mesh<Context>>
  {
    /// The execution context type
    using ContextType = Context;
  };
}

namespace Rodin::Geometry
{
  /**
   * @brief Mesh data structures, geometric operations, and computational geometry utilities.
   *
   * The Geometry module provides comprehensive support for mesh generation, manipulation,
   * and geometric computations in finite element analysis. It includes data structures
   * for representing meshes, connectivity information, geometric transformations, and
   * algorithms for mesh processing.
   *
   * # Key Components
   *
   * ## Mesh Data Structures
   * - **Mesh**: Core unstructured mesh class with support for mixed element types
   * - **SubMesh**: Mesh representing a subregion of a parent mesh
   * - **Shard**: Mesh shard for distributed parallel computing
   *
   * ## Connectivity Management
   * - Automatic computation of incidence relations @f$ d \rightarrow d' @f$
   * - Efficient storage using sparse data structures
   * - Support for arbitrary polytope dimensions
   *
   * ## Geometric Transformations
   * - Reference-to-physical element mappings
   * - Isoparametric transformations
   * - Jacobian computations
   *
   * ## Mesh Generation
   * - Structured mesh generation (uniform grids)
   * - Support for various file formats (MFEM, GMSH, MEDIT)
   * - Mesh refinement and coarsening
   *
   * ## Partitioning
   * - Domain decomposition for parallel computing
   * - Ghost polytope tracking
   * - Load balancing strategies
   *
   * # Mathematical Foundation
   *
   * A mesh @f$ \mathcal{T}_h @f$ is a partition of a domain @f$ \Omega @f$:
   * @f[
   *   \Omega = \bigcup_{K \in \mathcal{T}_h} K
   * @f]
   * where each element @f$ K @f$ is a polytope (simplex or hypercube).
   *
   * @see Mesh, Polytope, Connectivity
   */

  /**
   * @brief Connected component labeling result.
   *
   * Represents the result of connected component analysis on a mesh,
   * storing disjoint connected regions.
   */
  class CCL
  {
    public:
      /**
       * @brief Type representing a connected component.
       *
       * A component is a set of polytope indices that are connected.
       */
      using Component = FlatSet<Index>;

      /**
       * @brief Constructs CCL from a deque of components.
       * @param[in] dq Deque of connected components
       */
      CCL(std::deque<Component>&& dq)
        : m_components(std::move(dq))
      {}

      /**
       * @brief Gets all connected components.
       * @returns Deque of connected components
       */
      const std::deque<Component>& getComponents() const
      {
        return m_components;
      }

      /**
       * @brief Gets iterator to first component.
       */
      auto begin()
      {
        return m_components.begin();
      }

      /**
       * @brief Gets iterator past last component.
       */
      auto end()
      {
        return m_components.end();
      }

      /**
       * @brief Gets const iterator to first component.
       */
      auto begin() const
      {
        return m_components.begin();
      }

      /**
       * @brief Gets const iterator past last component.
       */
      auto end() const
      {
        return m_components.end();
      }

      /**
       * @brief Gets const iterator to first component.
       */
      auto cbegin() const
      {
        return m_components.cbegin();
      }

      /**
       * @brief Gets const iterator past last component.
       */
      auto cend() const
      {
        return m_components.cend();
      }

      /**
       * @brief Gets the number of connected components.
       * @returns Component count
       */
      size_t getCount() const
      {
        return m_components.size();
      }

    private:
      std::deque<Component> m_components;
  };

  /**
   * @defgroup MeshTypes Mesh Types and Template Specializations
   * @brief Different types of mesh and template specializations of the Mesh
   * class.
   * @see Mesh
   */

  /**
   * @ingroup RodinGeometry
   * @brief Abstract base class for all mesh implementations.
   *
   * MeshBase provides the common interface for all mesh types in Rodin,
   * defining fundamental operations such as:
   * - Point inclusion testing
   * - Basic mesh properties (dimension, vertex count, etc.)
   * - Polytope access through iterators
   * - Mesh I/O operations
   * - Connectivity queries
   *
   * # Inheritance Hierarchy
   *
   * All concrete mesh implementations derive from this base class:
   * - Mesh<Context::Local>: Sequential (non-distributed) meshes
   * - SubMesh<Context::Local>: Submeshes representing subregions
   * - Shard: Distributed mesh shards for parallel computing
   *
   * # Thread Safety
   * MeshBase and derived classes are not thread-safe during construction.
   * Once finalized, read-only operations are thread-safe.
   *
   * @see Mesh, SubMesh, Shard
   */
  class MeshBase
  {
    public:
      /**
       * @brief Virtual destructor.
       */
      virtual ~MeshBase() = default;

      /**
       * @brief Equality comparison (pointer equality).
       * @param[in] other Mesh to compare with
       * @returns True if this is the same mesh object
       */
      constexpr
      bool operator==(const MeshBase& other) const
      {
        return this == &other;
      }

      /**
       * @brief Inequality comparison.
       * @param[in] other Mesh to compare with
       * @returns True if not the same mesh object
       */
      constexpr
      bool operator!=(const MeshBase& other) const
      {
        return this != &other;
      }

      /**
       * @brief Tests point inclusion in the mesh.
       * @param[in] p Point to test
       * @returns Optional point if included in mesh
       *
       * Determines whether a point lies within the mesh domain and
       * returns its representation if so.
       */
      virtual Optional<Point> inclusion(const Point& p) const;

      /**
       * @brief Checks if mesh is empty.
       * @returns True if mesh has no vertices
       *
       * An empty mesh is defined as a mesh with no vertices.
       */
      bool isEmpty() const
      {
        return getVertexCount() == 0;
      }

      /**
       * @brief Checks if mesh is a surface mesh.
       * @returns True if mesh has codimension 1
       *
       * A mesh is considered a surface mesh if its codimension is 1,
       * meaning the difference between space dimension and topological
       * dimension equals 1.
       *
       * Examples:
       * - 2D mesh embedded in 3D space (triangulated surface)
       * - 1D mesh embedded in 2D space (curve)
       */
      bool isSurface() const;

      /**
       * @brief Gets the number of vertices.
       * @returns Vertex count
       */
      size_t getVertexCount() const
      {
        return getPolytopeCount(0);
      }

      /**
       * @brief Gets the number of faces.
       * @returns Face count (codimension-1 polytopes)
       */
      size_t getFaceCount() const
      {
        return getPolytopeCount(getDimension() - 1);
      }

      /**
       * @brief Gets the number of cells.
       * @returns Cell count (maximal dimension polytopes)
       */
      size_t getCellCount() const
      {
        return getPolytopeCount(getDimension());
      }

      /**
       * @brief Gets face attribute.
       * @param[in] index Face index
       * @returns Attribute value
       */
      Attribute getFaceAttribute(Index index) const
      {
        return getAttribute(getDimension() - 1, index);
      }

      /**
       * @brief Gets cell attribute.
       * @param[in] index Cell index
       * @returns Attribute value
       */
      Attribute getCellAttribute(Index index) const
      {
        return getAttribute(getDimension(), index);
      }

      /**
       * @brief Scales mesh coordinates by a factor.
       * @param[in] c Scaling factor
       * @returns Reference to this mesh
       *
       * Multiplies all vertex coordinates by @f$ c @f$.
       */
      virtual MeshBase& scale(Real c) = 0;

      /**
       * @brief Loads mesh from file.
       * @param[in] filename Path to mesh file
       * @param[in] fmt File format (MFEM, GMSH, etc.)
       * @returns Reference to this mesh
       */
      virtual MeshBase& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) = 0;

      /**
       * @brief Saves mesh to file.
       * @param[in] filename Path to output file
       * @param[in] fmt File format
       */
      virtual void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) const = 0;

      /**
       * @brief Flushes cached mesh data.
       *
       * Clears any internally cached derived data (e.g., cached Jacobians).
       */
      virtual void flush() = 0;

      /**
       * @brief Checks if this mesh is a submesh.
       * @returns True if this is a SubMesh instance
       *
       * A Mesh which is also a SubMesh may be casted down to access
       * SubMesh-specific functionality:
       * @code{.cpp}
       * if (mesh.isSubMesh())
       * {
       *   auto& submesh = static_cast<SubMesh<Context::Local>&>(mesh);
       *   const auto& parent = submesh.getParent();
       * }
       * @endcode
       */
      virtual bool isSubMesh() const = 0;

      /**
       * @brief Checks if a face is an interface between regions.
       * @param[in] faceIdx Face index
       * @returns True if face connects cells with different attributes
       */
      virtual bool isInterface(Index faceIdx) const = 0;

      /**
       * @brief Checks if a face is on the mesh boundary.
       * @param[in] faceIdx Face index
       * @returns True if face belongs to only one cell
       */
      virtual bool isBoundary(Index faceIdx) const = 0;

      /**
       * @brief Casts to SubMeshBase interface (mutable).
       * @returns Reference to SubMeshBase interface
       * @throws Exception if mesh is not a submesh
       */
      virtual SubMeshBase& asSubMesh() = 0;

      /**
       * @brief Casts to SubMeshBase interface (const).
       * @returns Const reference to SubMeshBase interface
       * @throws Exception if mesh is not a submesh
       */
      virtual const SubMeshBase& asSubMesh() const = 0;

      /**
       * @brief Gets the topological dimension.
       * @returns Mesh dimension @f$ D @f$ (highest polytope dimension)
       * @see getSpaceDimension()
       */
      virtual size_t getDimension() const = 0;

      /**
       * @brief Gets the ambient space dimension.
       * @returns Space dimension in which mesh is embedded
       * @see getDimension()
       *
       * For example, a 2D triangulated surface in 3D space has:
       * - getDimension() = 2
       * - getSpaceDimension() = 3
       */
      virtual size_t getSpaceDimension() const = 0;

      /**
       * @brief Gets the total mesh volume.
       * @returns Sum of all cell volumes
       */
      virtual Real getVolume() const = 0;

      /**
       * @brief Gets volume of cells with specified attribute.
       * @param[in] attr Attribute value
       * @returns Sum of volumes of cells with attribute @p attr
       * @note Returns 0 if attribute does not exist
       */
      virtual Real getVolume(Attribute attr) const = 0;

      /**
       * @brief Gets volume of cells with attributes in set.
       * @param[in] attr Set of attribute values
       * @returns Sum of volumes
       */
      virtual Real getVolume(const FlatSet<Attribute>& attr) const = 0;

      /**
       * @brief Gets the total mesh perimeter.
       * @returns Sum of all cell perimeters
       */
      virtual Real getPerimeter() const = 0;

      /**
       * @brief Gets perimeter of cells with specified attribute.
       * @param[in] attr Attribute value
       * @returns Sum of perimeters
       * @note Returns 0 if attribute does not exist
       */
      virtual Real getPerimeter(Attribute attr) const = 0;

      /**
       * @brief Gets perimeter of cells with attributes in set.
       * @param[in] attr Set of attribute values
       * @returns Sum of perimeters
       */
      virtual Real getPerimeter(const FlatSet<Attribute>& attr) const = 0;

      /**
       * @brief Gets total area.
       * @returns Sum of all element areas
       */
      virtual Real getArea() const = 0;

      /**
       * @brief Gets area of elements with specified attribute.
       * @param[in] attr Attribute value
       * @returns Sum of areas
       */
      virtual Real getArea(Attribute attr) const = 0;

      /**
       * @brief Gets area of elements with attributes in set.
       * @param[in] attr Set of attribute values
       * @returns Sum of areas
       */
      virtual Real getArea(const FlatSet<Attribute>& attr) const = 0;

      /**
       * @brief Gets measure of polytopes of given dimension.
       * @param[in] d Polytope dimension
       * @returns Total measure
       */
      virtual Real getMeasure(size_t d) const = 0;

      /**
       * @brief Gets measure of polytopes with attribute.
       * @param[in] d Polytope dimension
       * @param[in] attr Attribute value
       * @returns Total measure
       */
      virtual Real getMeasure(size_t d, Attribute attr) const = 0;

      /**
       * @brief Gets measure of polytopes with attributes in set.
       * @param[in] d Polytope dimension
       * @param[in] attr Set of attribute values
       * @returns Total measure
       */
      virtual Real getMeasure(size_t d, const FlatSet<Attribute>& attr) const = 0;

      /**
       * @brief Gets iterator over boundary faces.
       * @returns FaceIterator for boundary faces
       */
      virtual FaceIterator getBoundary() const = 0;

      /**
       * @brief Gets iterator over interface faces.
       * @returns FaceIterator for interface faces
       */
      virtual FaceIterator getInterface() const = 0;

      /**
       * @brief Gets count of polytopes of given dimension.
       * @param[in] dimension Polytope dimension
       * @returns Number of polytopes
       */
      virtual size_t getPolytopeCount(size_t dimension) const = 0;

      /**
       * @brief Gets count of polytopes of given geometry type.
       * @param[in] g Polytope geometry type
       * @returns Number of polytopes
       */
      virtual size_t getPolytopeCount(Polytope::Type g) const = 0;

      /**
       * @brief Gets iterator over all cells.
       * @returns CellIterator
       */
      virtual CellIterator getCell() const = 0;

      /**
       * @brief Gets iterator over all faces.
       * @returns FaceIterator
       */
      virtual FaceIterator getFace() const = 0;

      /**
       * @brief Gets iterator over all vertices.
       * @returns VertexIterator
       */
      virtual VertexIterator getVertex() const = 0;

      /**
       * @brief Gets iterator over polytopes of given dimension.
       * @param[in] dimension Polytope dimension
       * @returns PolytopeIterator
       */
      virtual PolytopeIterator getPolytope(size_t dimension) const = 0;

      /**
       * @brief Gets iterator to specific cell.
       * @param[in] idx Cell index
       * @returns CellIterator positioned at cell
       */
      virtual CellIterator getCell(Index idx) const = 0;

      /**
       * @brief Gets iterator to specific face.
       * @param[in] idx Face index
       * @returns FaceIterator positioned at face
       */
      virtual FaceIterator getFace(Index idx) const = 0;

      /**
       * @brief Gets iterator to specific vertex.
       * @param[in] idx Vertex index
       * @returns VertexIterator positioned at vertex
       */
      virtual VertexIterator getVertex(Index idx) const = 0;

      /**
       * @brief Gets iterator to specific polytope.
       * @param[in] dimension Polytope dimension
       * @param[in] idx Polytope index
       * @returns PolytopeIterator positioned at polytope
       */
      virtual PolytopeIterator getPolytope(size_t dimension, Index idx) const = 0;

      /**
       * @brief Gets polytope transformation.
       * @param[in] dimension Polytope dimension
       * @param[in] idx Polytope index
       * @returns Reference to transformation @f$ x: K \rightarrow \tau @f$
       */
      virtual const PolytopeTransformation& getPolytopeTransformation(size_t dimension, Index idx) const = 0;

      /**
       * @brief Gets geometry type of a polytope.
       * @param[in] dimension Polytope dimension
       * @param[in] idx Polytope index
       * @returns Geometry type
       */
      virtual Polytope::Type getGeometry(size_t dimension, Index idx) const = 0;

      /**
       * @brief Gets attribute of a polytope.
       * @param[in] dimension Polytope dimension
       * @param[in] index Polytope index
       * @returns Attribute value
       */
      virtual Attribute getAttribute(size_t dimension, Index index) const = 0;

      /**
       * @brief Sets attribute of a polytope.
       * @param[in] p Pair of (dimension, index)
       * @param[in] attr New attribute value
       * @returns Reference to this mesh
       */
      virtual MeshBase& setAttribute(const std::pair<size_t, Index>& p, Attribute attr) = 0;

      /**
       * @brief Gets reference to mesh connectivity.
       * @returns Const reference to Connectivity object
       */
      virtual ConnectivityBase& getConnectivity() = 0;

      /**
       * @brief Gets a constant reference to the mesh connectivity.
       */
      virtual const ConnectivityBase& getConnectivity() const = 0;

      /**
       * @brief Gets the space coordinates of the vertex at the given index.
       * @param[in] idx Vertex index
       */
      virtual Eigen::Map<const Math::SpatialPoint> getVertexCoordinates(Index idx) const = 0;

      /**
       * @brief Sets the space coordinate of the vertex at the given index for
       * the given coordinate index.
       * @param[in] idx Vertex index
       * @param[in] s New coordinate
       * @param[in] i Coordinate index
       *
       * For example, the following code sets the coordinates of the 0-vertex
       * to @f$ (0, 5, 10) @f$ in a mesh embedded in three dimensional space.
       *
       * @code{.cpp}
       *   Mesh mesh;
       *   // Add vertices...
       *   mesh.setVertexCoordinates(0, 0.0, 0);
       *   mesh.setVertexCoordinates(0, 5.0, 1);
       *   mesh.setVertexCoordinates(0, 10.0, 2);
       * @endcode
       */
      virtual MeshBase& setVertexCoordinates(Index idx, Real s, size_t i) = 0;

      /**
       * @brief Sets the space coordinate of the vertex at the given index for
       * the given coordinate index.
       * @param[in] idx Vertex index
       * @param[in] coords New coordinates
       */
      virtual MeshBase& setVertexCoordinates(Index idx, const Math::SpatialPoint& coords) = 0;

      virtual MeshBase& setPolytopeTransformation(
          const std::pair<size_t, Index> p, PolytopeTransformation* trans) = 0;

      /**
       * @brief Gets the execution context of the mesh.
       * @returns Reference to the context object
       *
       * The context determines whether the mesh is local, distributed, or
       * uses another execution model.
       */
      virtual const Context::Base& getContext() const = 0;
  };

  /// Type alias for Mesh<Context::Local>
  using LocalMesh =
    Mesh<Context::Local>;

  /**
   *
   * @ingroup MeshTypes
   * @brief Represents the subdivision of some domain into faces of (possibly)
   * different geometries.
   */
  template <>
  class Mesh<Context::Local> : public MeshBase
  {
    friend class boost::serialization::access;

    public:
      using Context =
        Context::Local;

      using Parent =
        MeshBase;

      /**
       * @brief Class used to build Mesh<Context::Local> instances.
       */
      class Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder()
            : m_initialized(false)
          {}

          virtual ~Builder() = default;

          /**
           * @brief Deleted copy constructor.
           */
          Builder(const Builder&) = delete;

          /**
           * @brief Move constructor.
           */
          Builder(Builder&&) = default;

          /**
           * @brief Deleted copy assignment.
           */
          Builder& operator=(const Builder&) = delete;

          /**
           * @brief Move assignment.
           */
          Builder& operator=(Builder&&) = default;

          /**
           * @brief Reserves memory accross the data structure for the polytopes of
           * the given dimension.
           */
          Builder& reserve(size_t d, size_t count);

          /**
           * @brief Initializes construction of the mesh object.
           */
          Builder& initialize(size_t sdim);

          /**
           * @brief Sets the number of nodes in the mesh.
           *
           * @note This method requires initialize(size_t) to be called
           * beforehand.
           */
          Builder& nodes(size_t n);

          /**
           * @brief Adds vertex with coordinates given by the fixed size array.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          template <size_t Size>
          Builder& vertex(const Real (&data)[Size])
          {
            assert(Size == m_sdim);
            m_vertices.col(m_nodes++) = data;
            return *this;
          }

          /**
           * @brief Adds vertex with coordinates given by the initializer list.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(std::initializer_list<Real> l);

          /**
           * @brief Adds vertex with coordinates given by the array pointer.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Real* data);

          /**
           * @brief Adds vertex with coordinates given by the mapped memory.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(Eigen::Map<const Math::SpatialPoint> x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(Math::Vector<Real>&& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Math::Vector<Real>& x);

          /**
           * @brief Sets the attribute of the given polytope.
           */
          Builder& attribute(const std::pair<size_t, Index>& p, Attribute attr);

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Type t, std::initializer_list<Index> vs)
          {
            return polytope(t, IndexArray({ vs }));
          }

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Type t, const IndexArray& vs);

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Type t, IndexArray&& vs);

          template <class T>
          Builder& segment(T&& vs)
          {
            return polytope(Polytope::Type::Segment, std::forward<T>(vs));
          }

          template <class T>
          Builder& quadrilateral(T&& vs)
          {
            return polytope(Polytope::Type::Quadrilateral, std::forward<T>(vs));
          }

          template <class T>
          Builder& triangle(T&& vs)
          {
            return polytope(Polytope::Type::Triangle, std::forward<T>(vs));
          }

          template <class T>
          Builder& tetrahedron(T&& vs)
          {
            return polytope(Polytope::Type::Tetrahedron, std::forward<T>(vs));
          }

          /**
           * @brief Finalizes construction of the Mesh<Context::Local> object.
           */
          Mesh finalize();

          /**
           * @brief Sets the vertex coordinates (copy).
           * @param[in] vertices Point matrix containing vertex coordinates
           * @returns Reference to this builder
           */
          Builder& setVertices(const Math::PointMatrix& vertices);

          /**
           * @brief Sets the vertex coordinates (move).
           * @param[in] vertices Point matrix containing vertex coordinates
           * @returns Reference to this builder
           */
          Builder& setVertices(Math::PointMatrix&& vertices);

          /**
           * @brief Sets the mesh connectivity.
           * @param[in] connectivity Connectivity object to move into the mesh
           * @returns Reference to this builder
           */
          Builder& setConnectivity(Connectivity<Context>&& connectivity);

          /**
           * @brief Sets the attribute index.
           * @param[in] attrIndex Attribute index object to move into the mesh
           * @returns Reference to this builder
           */
          Builder& setAttributeIndex(AttributeIndex&& attrIndex);

          /**
           * @brief Sets the transformation index.
           * @param[in] transIndex Transformation index object to move into the mesh
           * @returns Reference to this builder
           */
          Builder& setTransformationIndex(PolytopeTransformationIndex&& transIndex);

          Connectivity<Context>& getConnectivity()
          {
            return m_connectivity;
          }

          const Connectivity<Context>& getConnectivity() const
          {
            return m_connectivity;
          }

        private:
          bool m_initialized;

          size_t m_sdim;
          size_t m_nodes;

          Math::PointMatrix m_vertices;
          Connectivity<Context> m_connectivity;

          AttributeIndex m_attributes;
          PolytopeTransformationIndex m_transformations;
      };

      /**
       * @brief Generates a Builder instance to build a Mesh object.
       */
      static Builder Build()
      {
        return Builder();
      }

      static Mesh UniformGrid(Polytope::Type g, std::initializer_list<size_t> l)
      {
        Array<size_t> shape(l.size());
        std::copy(l.begin(), l.end(), shape.begin());
        return UniformGrid(g, shape);
      }

      static Mesh Box(Polytope::Type g, std::initializer_list<size_t> l)
      {
        Array<size_t> shape(l.size());
        std::copy(l.begin(), l.end(), shape.begin());
        return Box(g, shape);
      }

      /**
       * @brief Generates a uniform grid for a given geometry.
       *
       * @param[in] g Geometry type of the cells (e.g., Triangle,
       * Quadrilateral)
       * @param[in] shape Number of nodes in each coordinate direction
       */
      static Mesh UniformGrid(Polytope::Type g, const Array<size_t>& shape);

      static Mesh Box(Polytope::Type g, const Array<size_t>& shape);

      /**
      * @brief Constructs an empty mesh with no cells.
      */
      Mesh()
        : m_sdim(0)
      {}

      Mesh(const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM)
      {
        load(filename, fmt);
      }

      /**
      * @brief Performs a copy of another mesh.
      */
      Mesh(const Mesh& other);

      /**
      * @brief Move constructs the mesh from another mesh.
      */
      Mesh(Mesh&& other);

      virtual ~Mesh() = default;

      Mesh& operator=(const Mesh& other) = delete;

      /**
      * @brief Move assigns the mesh from another mesh.
      */
      Mesh& operator=(Mesh&&);

      /**
       * @brief Displaces the mesh nodes by the displacement @f$ u @f$.
       * @param[in] u Displacement at each node
       *
       * Given a vector valued function @f$ \vec{u} @f$, the method will perform the
       * displacement
       * @f[
       *   x \mapsto x + \vec{u}(x)
       * @f]
       * at each node @f$ x @f$ of the mesh.
       *
       * @note The vector dimension of @f$ \vec{u} @f$ must be equal to the
       * space dimension.
       *
       * @returns Reference to this (for method chaining)
       */
      template <class FunctionDerived>
      Mesh& displace(const Variational::FunctionBase<FunctionDerived>& u)
      {
        for (auto it = getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(
              *it,
              Polytope::Traits(Polytope::Type::Point).getVertex(0),
              it->getCoordinates());
          m_vertices.col(it->getIndex()) += u(p);
        }
        return *this;
      }


      virtual void flush() override
      {
        m_transformations.clear();
      }

      /**
       * @brief Gets the total volume of the mesh.
       * @returns Sum of all cell volumes.
       */
      Real getVolume() const override;

      /**
       * @brief Gets the sum of the volumes of the cells given by the
       * specified attribute.
       * @param[in] attr Attribute of cells
       * @returns Sum of element volumes with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the volume.
       */
      Real getVolume(Attribute attr) const override;

      Real getVolume(const FlatSet<Attribute>& attr) const override;

      /**
       * @brief Gets the total perimeter of the mesh.
       * @returns Sum of all element perimeters.
       */
      Real getPerimeter() const override;

      /**
       * @brief Gets the sum of the perimeters of the cells given by the
       * specified attribute.
       * @param[in] attr Attribute of cells
       * @returns Sum of element perimeters with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the perimeter.
       */
      Real getPerimeter(Attribute attr) const override;

      Real getPerimeter(const FlatSet<Attribute>& attr) const override;

      Real getArea() const override;

      Real getArea(Attribute attr) const override;

      Real getArea(const FlatSet<Attribute>& attr) const override;

      Real getMeasure(size_t d) const override;

      Real getMeasure(size_t d, Attribute attr) const override;

      Real getMeasure(size_t d, const FlatSet<Attribute>& attr) const override;

      template <class BinaryPredicate>
      CCL ccl(const BinaryPredicate& p) const
      {
        return ccl(getDimension(), p, [](const Polytope&) { return true; });
      }

      template <class BinaryPredicate>
      CCL ccl(size_t d, const BinaryPredicate& p) const
      {
        return ccl(d, p, [](const Polytope&) { return true; });
      }

      template <class BinaryPredicate>
      CCL ccl(size_t d, const BinaryPredicate& p, Attribute attr) const
      {
        return ccl(d, p,
          [attr](const Polytope& polytope) { return attr == polytope.getAttribute(); });
      }

      template <class BinaryPredicate>
      CCL ccl(size_t d, const BinaryPredicate& p, const FlatSet<Attribute>& attrs) const
      {
        return ccl(d, p,
          [&attrs](const Polytope& polytope) { return attrs.size() == 0 || attrs.contains(polytope.getAttribute()); });
      }

      template <class BinaryPredicate, class UnitaryPredicate>
      CCL ccl(size_t d, const BinaryPredicate& p, const UnitaryPredicate& f) const
      {
        FlatSet<Index> visited;
        visited.reserve(getPolytopeCount(d));
        std::deque<Index> searchQueue;
        std::deque<FlatSet<Index>> res;

        // Perform the labelling
        for (auto it = getPolytope(d); it; ++it)
        {
          const Index i = it->getIndex();
          if (!visited.count(i))
          {
            if (f(*it))
            {
              res.push_back({});
              searchQueue.push_back(i);
            }
            while (searchQueue.size() > 0)
            {
              const Index idx = searchQueue.back();
              const auto el = getPolytope(d, idx);
              searchQueue.pop_back();
              const auto result = visited.insert(idx);
              const Boolean inserted = result.second;
              if (inserted)
              {
                res.back().insert(idx);
                for (auto adj = el->getAdjacent(); adj; ++adj)
                {
                  if (p(*el, *adj))
                  {
                    if (f(*adj))
                      searchQueue.push_back(adj->getIndex());
                  }
                }
              }
            }
          }
        }
        return res;
      }

      /**
      * @brief Skins the mesh to obtain its boundary mesh
      * @returns SubMesh object to the boundary region of the mesh
      *
      * This function "skins" the mesh to return the boundary as a new SubMesh
      * object. The resulting mesh will be embedded in the original space
      * dimension.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 2
      * @f$ are included if the connectivity @f$ (D - 1) \longrightarrow d @f$
      * is already computed in the mesh.
      */
      virtual SubMesh<Context> skin() const;

      /**
      * @brief Trims the cells with the given attribute.
      * @param[in] attr Attribute to trim
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call trim(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context> trim(Attribute attr) const;

      /**
      * @brief Trims the cells with the given attribute.
      * @param[in] attrs Attributes to trim
      * @returns SubMesh object to the remaining region of the mesh
      *
      * This function will trim discard all the cells that have an attribute
      * in the given set of attributes.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 1
      * @f$ are included if the connectivity @f$ D \longrightarrow d @f$ is
      * already computed in the mesh.
      *
      * @returns A SubMesh object consisting of cells that have attributes
      * not in the given set.
      */
      virtual SubMesh<Context> trim(const FlatSet<Attribute>& attrs) const;

      /**
      * @brief Keeps the cells with the given attribute.
      * @param[in] attr Attribute to keep
      *
      * Convenience function to call keep(const std::FlatSet<Attribute>&) with
      * only one attribute.
      *
      * @returns SubMesh of the remaining region mesh
      */
      virtual SubMesh<Context> keep(Attribute attr) const;

      /**
      * @brief Trims the cells with the given attributes.
      * @param[in] attrs Attributes to trim
      *
      * This function will trim keep only the cells that have an attribute
      * in the given set of attributes.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 1
      * @f$ are included if the connectivity @f$ D \longrightarrow d @f$ is
      * already computed in the mesh.
      *
      * @returns A SubMesh object consisting of cells that have attributes
      * not in the given set.
      */
      virtual SubMesh<Context> keep(const FlatSet<Attribute>& attrs) const;

      Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap)
      {
        return trace(tmap, FlatSet<Attribute>{});
      }

      Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap, Attribute attr)
      {
        return trace(tmap, FlatSet<Attribute>{ attr });
      }

      virtual Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap, const FlatSet<Attribute>& attrs);

      Mesh& trace(const Map<Attribute, Attribute>& tmap)
      {
        return trace(tmap, FlatSet<Attribute>{});
      }

      virtual Mesh& trace(const Map<Attribute, Attribute>& tmap, const FlatSet<Attribute>& attrs);

      SubMeshBase& asSubMesh() override;

      const SubMeshBase& asSubMesh() const override;

      /**
      * @brief Loads a mesh from file in the given format.
      * @param[in] filename Name of file to read
      * @param[in] fmt Mesh file format
      * @returns Reference to this (for method chaining)
      */
      virtual Mesh& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) override;

      /**
      * @brief Saves a mesh to file in the given format.
      * @param[in] filename Name of file to write
      * @param[in] fmt Mesh file format
      */
      virtual void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) const override;

      virtual Mesh& scale(Real c) override;

      const AttributeIndex& getAttributeIndex() const
      {
        return m_attributes;
      }

      const PolytopeTransformationIndex& getPolytopeTransformationIndex() const
      {
        return m_transformations;
      }

      const Math::PointMatrix& getVertices() const
      {
        return m_vertices;
      }

      const Context& getContext() const override
      {
        return m_context;
      }

      virtual size_t getPolytopeCount(size_t dim) const override;

      virtual size_t getPolytopeCount(Polytope::Type g) const override;

      virtual CellIterator getCell() const override;

      virtual FaceIterator getFace() const override;

      virtual VertexIterator getVertex() const override;

      virtual PolytopeIterator getPolytope(size_t dimension) const override;

      virtual FaceIterator getBoundary() const override;

      virtual FaceIterator getInterface() const override;

      virtual CellIterator getCell(Index idx) const override;

      virtual FaceIterator getFace(Index idx) const override;

      virtual VertexIterator getVertex(Index idx) const override;

      virtual PolytopeIterator getPolytope(size_t dimension, Index idx) const override;

      virtual bool isSubMesh() const override
      {
        return false;
      }

      virtual bool isInterface(Index faceIdx) const override;

      virtual bool isBoundary(Index faceIdx) const override;

      virtual size_t getDimension() const override;

      virtual size_t getSpaceDimension() const override;

      virtual Polytope::Type getGeometry(size_t dimension, Index idx) const override;

      virtual Attribute getAttribute(size_t dimension, Index index) const override;

      virtual Connectivity<Context>& getConnectivity() override
      {
        return m_connectivity;
      }

      virtual const Connectivity<Context>& getConnectivity() const override
      {
        return m_connectivity;
      }

      virtual Eigen::Map<const Math::SpatialPoint> getVertexCoordinates(Index idx) const override;

      virtual Mesh& setAttribute(const std::pair<size_t, Index>&, Attribute attr) override;

      virtual Mesh& setVertexCoordinates(Index idx, Real xi, size_t i) override;

      virtual Mesh& setVertexCoordinates(Index idx, const Math::SpatialPoint& coords) override;

      virtual Mesh& setPolytopeTransformation(
          const std::pair<size_t, Index> p, PolytopeTransformation* trans) override;

      virtual const PolytopeTransformation& getPolytopeTransformation(
          size_t dimension, Index idx) const override;

      virtual PolytopeTransformation* getDefaultPolytopeTransformation(size_t d, Index i) const;

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & m_sdim;
        ar & m_vertices;
        ar & m_connectivity;
        ar & m_transformations;
        ar & m_attributes;
        ar & m_context;
      }

    private:
      size_t m_sdim;

      Math::PointMatrix m_vertices;
      Connectivity<Context> m_connectivity;

      AttributeIndex m_attributes;
      PolytopeTransformationIndex m_transformations;

      Context m_context;
  };
}

#endif
