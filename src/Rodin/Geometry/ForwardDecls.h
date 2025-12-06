/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_FORWARDDECLS_H
#define RODIN_MESH_FORWARDDECLS_H

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
#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Template container indexed by polytope geometry types.
   * @tparam T Type of values to store
   * @see GeometryIndexed.h
   */
  template <class T>
  class GeometryIndexed;

  /**
   * @brief Index generator that wraps a container's iterators.
   * @tparam T Container type
   * @see IndexGenerator.h
   */
  template <class T>
  class ContainerIndexGenerator;

  /**
   * @brief Enumeration of polytope geometries (forward declaration).
   * @see Polytope::Type
   */
  enum class Type;

  /**
   * @brief Base class for polytope transformations.
   * @see PolytopeTransformation.h
   */
  class PolytopeTransformation;

  /**
   * @brief Template for isoparametric transformations.
   * @tparam FE Finite element type
   * @see IsoparametricTransformation.h
   */
  template <class FE>
  class IsoparametricTransformation;

  /**
   * @brief Represents a geometric polytope in a mesh.
   * @see Polytope.h
   */
  class Polytope;

  /**
   * @brief Represents a mesh cell (polytope of maximal dimension).
   * @see Polytope.h
   */
  class Cell;

  /**
   * @brief Represents a mesh face (polytope of dimension d-1).
   * @see Polytope.h
   */
  class Face;

  /**
   * @brief Represents a mesh vertex (0-dimensional polytope).
   * @see Polytope.h
   */
  class Vertex;

  /**
   * @brief Represents a point on a mesh polytope.
   * @see Point.h
   */
  class Point;

  /**
   * @brief Iterator over polytopes.
   * @see PolytopeIterator.h
   */
  class PolytopeIterator;

  /**
   * @brief Iterator over mesh cells.
   * @see PolytopeIterator.h
   */
  class CellIterator;

  /**
   * @brief Iterator over mesh faces.
   * @see PolytopeIterator.h
   */
  class FaceIterator;

  /**
   * @brief Iterator over mesh vertices.
   * @see PolytopeIterator.h
   */
  class VertexIterator;

  /**
   * @brief Base class for mesh representations.
   * @see Mesh.h
   */
  class MeshBase;

  /**
   * @brief Template for connectivity information.
   * @tparam ContextType Context type (Local, MPI, etc.)
   * @see Connectivity.h
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
   * The template parameter @p ContextType determines the execution model:
   * - Context::Local for sequential execution
   * - Context::MPI for distributed parallel execution
   *
   * @tparam ContextType Execution context (default: Context::Local)
   * @see MeshBase, Mesh.h
   */
  template <class ContextType = Context::Local>
  class Mesh;

  /**
   * @brief Base class for SubMesh functionality.
   * @see SubMesh.h
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
   * @tparam Context Execution context type
   * @see SubMeshBase, SubMesh.h
   */
  template <class Context>
  class SubMesh;

  /**
   * @brief Sequential SubMesh specialization.
   * @see SubMesh.h
   */
  template <>
  class SubMesh<Context::Local>;

  /**
   * @brief Builder for constructing SubMesh objects.
   * @tparam Context Execution context type
   * @see SubMesh.h
   */
  template <class Context>
  class SubMeshBuilder;
}

#endif
