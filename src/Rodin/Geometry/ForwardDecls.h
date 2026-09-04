/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_FORWARDDECLS_H
#define RODIN_GEOMETRY_FORWARDDECLS_H

/**
 * @file
 * @brief Forward declarations for the Rodin::Geometry module.
 *
 * This file provides forward declarations for key classes in the geometry
 * module to minimize compilation dependencies. Detailed documentation for
 * each class can be found in its respective header file.
 */

#include <cstddef>

#include "Rodin/Context.h"
#ifdef RODIN_USE_MPI
#include "Rodin/MPI/Context/ForwardDecls.h"
#endif
#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Template container indexed by polytope geometry types.
   * @tparam T Type of values to store
   * @see <a href="_geometry_indexed_8h.html">GeometryIndexed.h</a>
   */
  template <class T>
  class GeometryIndexed;

  /**
   * @brief Index generator that wraps a container's iterators.
   * @tparam T Container type
   * @see <a href="_index_generator_8h.html">IndexGenerator.h</a>
   */
  template <class T>
  class ContainerIndexGenerator;

  /**
   * @brief Enumeration of polytope geometries (forward declaration).
   * @see <a href="class_rodin_1_1_geometry_1_1_polytope.html#a1d1cfd8ffb84e947f82999c682b666a7">Polytope::Type</a>
   */
  enum class Type;

  /**
   * @brief Base class for polytope transformations.
   * @see <a href="_polytope_transformation_8h.html">PolytopeTransformation.h</a>
   */
  class PolytopeTransformation;

  /**
   * @brief Template for finite-element parametric transformations.
   * @tparam FE Finite element type
   * @see <a href="_parametric_transformation_8h.html">ParametricTransformation.h</a>
   */
  template <class FE>
  class ParametricTransformation;

  /**
   * @brief Represents a geometric polytope in a mesh.
   * @see <a href="_polytope_8h.html">Polytope.h</a>
   */
  class Polytope;

  class PolytopeQuadrature;

  /**
   * @brief Represents a mesh cell (polytope of maximal dimension).
   * @see <a href="_polytope_8h.html">Polytope.h</a>
   */
  class Cell;

  /**
   * @brief Represents a mesh face (polytope of dimension d-1).
   * @see <a href="_polytope_8h.html">Polytope.h</a>
   */
  class Face;

  /**
   * @brief Represents a mesh vertex (0-dimensional polytope).
   * @see <a href="_polytope_8h.html">Polytope.h</a>
   */
  class Vertex;

  /**
   * @brief Represents a point on a mesh polytope.
   * @see <a href="_point_8h.html">Point.h</a>
   */
  class Point;

  /**
   * @brief Iterator over polytopes.
   * @see <a href="_polytope_iterator_8h.html">PolytopeIterator.h</a>
   */
  class PolytopeIterator;

  /**
   * @brief Iterator over mesh cells.
   * @see <a href="_polytope_iterator_8h.html">PolytopeIterator.h</a>
   */
  class CellIterator;

  /**
   * @brief Iterator over mesh faces.
   * @see <a href="_polytope_iterator_8h.html">PolytopeIterator.h</a>
   */
  class FaceIterator;

  /**
   * @brief Iterator over mesh vertices.
   * @see <a href="_polytope_iterator_8h.html">PolytopeIterator.h</a>
   */
  class VertexIterator;

  /**
   * @brief Base class for mesh representations.
   * @see <a href="_geometry_2_mesh_8h.html">Mesh.h</a>
   */
  class MeshBase;

  class ConnectivityBase;

  /**
   * @brief Template for connectivity information.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Connectivity "Connectivity<Context::Local>" | Sequential connectivity table storing all requested incidence relations for a local mesh. |
   * | @ref Connectivity "Connectivity<Context::MPI>" | Distributed connectivity facade for rank-local shard topology and MPI mesh operations. |
   *
   * @tparam ContextType Context type (Local, MPI, etc.)
   * @see <a href="_connectivity_8h.html">Connectivity.h</a>
   */
  template <class ContextType>
  class Connectivity;

  /**
   * @brief Represents a polyhedral complex.
   *
   * A Mesh object represents a polyhedral complex @f$ \mathcal{T}_h @f$,
   * which is a set containing finitely many convex polyhedra. The mesh
   * provides access to polytopes of various dimensions and their connectivity.
   *
   * # Key Features
   * - Access to polytopes by dimension (vertices, faces, cells)
   * - Connectivity information between polytopes
   * - Support for attributes marking different regions
   * - Thread-safe operations (context-dependent)
   *
   * # Context Types
   * The template parameter @p ContextType determines the execution model.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | [Mesh<Context::Local>](class_rodin_1_1_geometry_1_1_mesh_3_01_context_1_1_local_01_4.html) | Sequential mesh storing the full incidence complex and geometry in a single process. |
   * | [Mesh<Context::MPI>](class_rodin_1_1_geometry_1_1_mesh_3_01_context_1_1_m_p_i_01_4.html) | Distributed mesh storing a rank-local shard, global/local index maps, and MPI context for parallel assembly. |
   *
   * @tparam ContextType Execution context (default: Context::Local)
   * @see <a href="class_rodin_1_1_geometry_1_1_mesh_base.html">MeshBase</a>
   * @see <a href="_geometry_2_mesh_8h.html">Mesh.h</a>
   */
  template <class ContextType = Context::Local>
  class Mesh;

  /**
   * @brief Base class for SubMesh functionality.
   * @see <a href="_geometry_2_sub_mesh_8h.html">SubMesh.h</a>
   */
  class SubMeshBase;

  /**
   * @brief Represents a subset of a Mesh.
   *
   * # Overview
   * A SubMesh is a mesh that represents a subset of another mesh, typically
   * corresponding to a specific region or boundary. It maintains references
   * to its parent mesh and provides mappings between polytope indices in
   * the child and parent meshes.
   *
   * # Mapping Between SubMesh and Parent Mesh
   * A SubMesh object holds a reference to its parent Mesh and includes
   * information about how polytopes and vertices are mapped between the
   * child and parent. This allows finite element operations on the SubMesh
   * to properly reference degrees of freedom in the parent mesh.
   *
   * # Downcasting
   * A Mesh that is also a SubMesh can be downcasted to access SubMesh
   * functionality. For instance:
   * @code{.cpp}
   * if (mesh.isSubMesh())
   * {
   *   // The cast is well defined
   *   auto& submesh = static_cast<SubMesh&>(mesh);
   *   const auto& parent = submesh.getParent();
   * }
   * @endcode
   *
   * # Use Cases
   * - Boundary condition application
   * - Interface problem formulation
   * - Region-specific operations
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | [SubMesh<Context::Local>](class_rodin_1_1_geometry_1_1_sub_mesh_3_01_context_1_1_local_01_4.html) | Sequential submesh view extracted from a local parent mesh. |
   * | [SubMesh<Context::MPI>](class_rodin_1_1_geometry_1_1_sub_mesh_3_01_context_1_1_m_p_i_01_4.html) | Distributed submesh view extracted from an MPI parent mesh and its rank-local shards. |
   *
   * @tparam Context Execution context type
   * @see <a href="class_rodin_1_1_geometry_1_1_sub_mesh_base.html">SubMeshBase</a>
   * @see <a href="_geometry_2_sub_mesh_8h.html">SubMesh.h</a>
   */
  template <class Context>
  class SubMesh;

  /**
   * @brief Sequential SubMesh specialization.
   * @see <a href="_geometry_2_sub_mesh_8h.html">SubMesh.h</a>
   */
  template <>
  class SubMesh<Context::Local>;

  /**
   * @brief Distributed MPI SubMesh specialization.
   * @see <a href="_m_p_i_2_geometry_2_sub_mesh_8h.html">Rodin/MPI/Geometry/SubMesh.h</a>
   */
#ifdef RODIN_USE_MPI
  template <>
  class SubMesh<Context::MPI>;
#endif

  /**
   * @brief Builder for constructing SubMesh objects.
   * @tparam Context Execution context type
   * @see <a href="_geometry_2_sub_mesh_8h.html">SubMesh.h</a>
   */
  template <class Context>
  class SubMeshBuilder;
}

#endif
