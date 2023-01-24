/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CONNECTIVITY_H
#define RODIN_GEOMETRY_CONNECTIVITY_H

#include <vector>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
   class Connectivity
   {
      public:
         Connectivity(size_t d, size_t dp)
            : m_d(d), m_dp(dp)
         {}

         void build(size_t d);

         void transpose(size_t d, size_t dp);

         void intersection(size_t d, size_t dp);

      private:
         size_t m_d, m_dp;
         std::vector<Index> m_indices;
         std::vector<size_t> m_offsets;
   };
}

#endif
