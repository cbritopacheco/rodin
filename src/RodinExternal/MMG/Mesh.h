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

#include "Configure.h"

namespace Rodin::External::MMG
{
   template <int Dimension, class Derived>
   class Mesh
   {
      public:
         virtual ~Mesh() = default;

         size_t getMaximumMemory() const
         {
            return getHandle()->memMax;
         }

         size_t getCurrentMemory() const
         {
            return getHandle()->memCur;
         }

         constexpr int getDimension() const
         {
            return Dimension;
         }

         MMG5_pMesh& getHandle()
         {
            return static_cast<Derived*>(this)->getHandle();
         }

         const MMG5_pMesh& getHandle() const
         {
            return static_cast<const Derived*>(this)->getHandle();
         }
   };
}

#endif
