/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_FORWARDDECLS_H
#define RODIN_MESH_FORWARDDECLS_H

#include <utility>

#include "Rodin/Traits.h"
#include "Rodin/Configure.h"

namespace Rodin
{
   namespace Context
   {
      struct Serial;
      struct Parallel;
   }
}

namespace Rodin::Geometry
{
   using Index = std::size_t;
   using Attribute = std::size_t;
   using Dimension = std::size_t;

   enum class Type;

   enum class Region
   {
      Domain,
      Boundary,
      Interface
   };

   class Simplex;

   class Element;

   class Face;

   class Vertex;

   class Interface;

   class Boundary;

   class Point;

   class SimplexVertexIterator;

   class AdjacentSimplexIterator;

   class AdjacentElementIterator;

   class AdjacentFaceIterator;

   class FaceElementIterator;

   class MeshSimplexIterator;

   class MeshElementIterator;

   class MeshFaceIterator;

   class MeshInterfaceIterator;

   class MeshBoundaryIterator;

   class MeshVertexIterator;

   class PointIterator;

   class MeshBase;

   /**
    * @brief Templated class for Mesh.
    * @tparam Trait Indicates whether if Mesh is in a parallel context. It is
    * one of Traits::Serial or Traits::Parallel.
    *
    * The Mesh class represents an n-dimensional support for instances of type
    * GridFunctionBase or ShapeFunctionBase.
    *
    * There are two possible specializations:
    * - Mesh<Traits::Serial>
    * - Mesh<Traits::Parallel>
    */
   template <class Trait = Context::Serial>
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
   template <class Trait>
   class SubMesh;
}

#endif
