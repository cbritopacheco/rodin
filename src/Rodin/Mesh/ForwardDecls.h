/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_FORWARDDECLS_H
#define RODIN_MESH_FORWARDDECLS_H

#include "Rodin/Traits.h"
#include "Rodin/Configure.h"

namespace Rodin
{
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
   template <class Trait = Traits::Serial>
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

   class ElementBase;

   class Element;
   class ElementView;

   class Face;
   class FaceView;
}

#endif
