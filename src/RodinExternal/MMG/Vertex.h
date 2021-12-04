/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_VERTEX_H
#define RODIN_RODININTEGRATION_MMG_VERTEX_H

#include <mmg/libmmg.h>

namespace Rodin::External::MMG
{
   template <int Dimension, class Derived>
   class Vertex
   {
      public:
         MMG5_pPoint& getHandle()
         {
            return Derived::getHandle();
         }

         const MMG5_pPoint& getHandle() const
         {
            return Derived::getHandle();
         }
   };
}

#endif
