/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_BASEMESH_H
#define RODIN_RODININTEGRATION_MMG_BASEMESH_H

#include <cassert>
#include <mmg/libmmg.h>

#include "Macros.h"

namespace Rodin::External::MMG
{
   template <class Derived>
   class BaseMesh
   {
      public:
         size_t getMaximumMemory() const
         {
            return getHandle()->memMax;
         }

         size_t getCurrentMemory() const
         {
            return getHandle()->memCur;
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
