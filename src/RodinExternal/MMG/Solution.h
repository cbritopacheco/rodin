/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_SOLUTION_H
#define RODIN_RODININTEGRATION_MMG_SOLUTION_H

#include <cassert>

#include <mmg/libmmg.h>
#include <mmg/mmg2d/libmmg2d.h>

namespace Rodin::External::MMG
{
   class BaseSolution
   {
      public:
         virtual ~BaseSolution() = default;
         virtual MMG5_pSol& getHandle() = 0;
         virtual const MMG5_pSol& getHandle() const = 0;
   };

   template <int Dimension, class T>
   class Solution : public BaseSolution
   {
      public:
         virtual ~Solution() = default;
         virtual MMG5_pSol& getHandle() = 0;
         virtual const MMG5_pSol& getHandle() const = 0;

         int getDimension() const
         {
            assert(getHandle()->dim == Dimension);
            return Dimension;
         }

         /**
          * @returns Number of points of the solution.
          */
         int count() const
         {
            return getHandle()->np;
         }

         /**
          * @returns Number of solutions per entity.
          */
         int size() const
         {
            return getHandle()->size;
         }
   };
}

#endif
