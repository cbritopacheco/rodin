/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_EDGE_H
#define RODIN_RODININTEGRATION_MMG_EDGE_H

#include "libmmg.h"

namespace Rodin::External::MMG
{
   class Edge
   {
      public:
         Edge(MMG5_pEdge edge);

      private:
         MMG5_pEdge m_edge;
   };
}

#endif
