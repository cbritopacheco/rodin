/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESH_HPP
#define RODIN_MESH_MESH_HPP

#include "Mesh.h"

#include "Element.h"

namespace Rodin::Geometry
{
   template <class T>
   std::enable_if_t<std::is_same_v<Element, T>, Element>
   MeshBase::get(int i) const
   {
      return Element(*this, getHandle().GetElement(i), i);
   }

   template <class T>
   std::enable_if_t<
      std::is_same_v<Element, T>,
      ElementView>
   MeshBase::get(int i)
   {
      return ElementView(*this, getHandle().GetElement(i), i);
   }

   template <class T>
   std::enable_if_t<std::is_same_v<BoundaryElement, T>,
      BoundaryElement>
   MeshBase::get(int i) const
   {
      return BoundaryElement(*this, getHandle().GetBdrElement(i), i);
   }

   template <class T>
   std::enable_if_t<
      std::is_same_v<BoundaryElement, T>,
      BoundaryElementView>
   MeshBase::get(int i)
   {
      return BoundaryElementView(*this, getHandle().GetBdrElement(i), i);
   }

   template <class T>
   std::enable_if_t<std::is_same_v<Face, T>,
      Face>
   MeshBase::get(int i) const
   {
      return Face(*this, getHandle().GetFace(i), i);
   }
}

#endif
