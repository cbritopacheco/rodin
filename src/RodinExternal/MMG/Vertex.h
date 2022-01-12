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
   template <int Dimension>
   class Vertex
   {
      public:
         virtual MMG5_pPoint& getHandle() = 0;
         virtual const MMG5_pPoint& getHandle() const = 0;
   };
}

#endif
