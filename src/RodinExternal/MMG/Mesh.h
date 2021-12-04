/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_MESH_H
#define RODIN_RODININTEGRATION_MMG_MESH_H

#include <cassert>
#include <mmg/libmmg.h>

#include "BaseMesh.h"

namespace Rodin::External::MMG
{
   template <int Dimension, class Derived>
   class Mesh : public BaseMesh<Mesh<Dimension, Derived>>
   {
      public:
         constexpr int getDimension() const
         {
            return Dimension;
         }

         MMG5_pMesh& getHandle()
         {
            return Derived::getHandle();
         }

         const MMG5_pMesh& getHandle() const
         {
            return Derived::getHandle();
         }
   };
}

#endif
