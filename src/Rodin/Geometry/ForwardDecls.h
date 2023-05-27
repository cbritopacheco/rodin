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
   * @brief Templated class for Mesh.
   */
  template <class ContextType = Context::Serial>
  class Mesh;

  /**
   * @brief Templated class for SubMesh.
   *
   * @tparam Trait Indicates whether the Mesh is in a parallel context. It is
   * one of Traits::Serial or Traits::Parallel.
   *
   * There are two possible specializations:
   * - SubMesh<Traits::Serial>
   * - SubMesh<Traits::Parallel>
   */
  template <class Context>
  class SubMesh;

  template <class Context>
  class SubMeshBuilder;
}

#endif
