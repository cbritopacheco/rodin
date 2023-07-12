/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_FORWARDDECLS_H
#define RODIN_MESH_FORWARDDECLS_H

#include <cstddef>

#include "Rodin/Context.h"
#include "Types.h"

namespace Rodin::Geometry
{
  template <class T>
  class ContainerIndexGenerator;

  enum class Type;

  class PolytopeTransformation;

  template <class FE>
  class IsoparametricTransformation;

  class Polytope;

  class Element;

  class Face;

  class Vertex;

  class Point;

  class PolytopeIterator;

  class ElementIterator;

  class FaceIterator;

  class VertexIterator;

  class MeshBase;

  /**
   * @brief Represents a polyhedral complex.
   *
   * A Mesh object represents a polyhedral complex @f$ \mathcal{T}_h @f$ is a
   * set containing finitely many convex polyhedra.
   */
  template <class ContextType = Context::Serial>
  class Mesh;

  /**
   * @brief Represents a subset of a Mesh.
   *
   * # Mapping between the SubMesh and the parent Mesh
   *
   * A SubMesh object holds a reference to its parent Mesh object and includes
   * details about how polytopes and vertices are mapped between the child and
   * parent Mesh.
   *
   * # Downcasting
   *
   * A Mesh that is also a SubMesh can be downcasted to access the SubMesh
   * functionality. For instance:
   * @code{.cpp}
   * if (mesh.isSubMesh())
   * {
   *   // The cast is well defined
   *   auto& submesh = static_cast<SubMesh&>(mesh);
   * }
   * @endcode
   */
  template <class Context>
  class SubMesh;

  template <class Context>
  class SubMeshBuilder;
}

#endif
