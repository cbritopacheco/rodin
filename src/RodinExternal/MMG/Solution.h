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
   /**
    * @brief Base class for MMG solution objects.
    * @tparam Dimension Mesh dimension
    *
    * This class wraps the functionality of the type MMG5_Sol.
    */
   class SolutionBase
   {
      public:
         virtual ~SolutionBase() = default;

         /**
          * @internal
          * @returns Reference to underlying solution handle.
          */
         virtual MMG5_pSol& getHandle() = 0;

         /**
          * @internal
          * @returns Constant reference to underlying solution handle.
          */
         virtual const MMG5_pSol& getHandle() const = 0;

         /**
          * @brief Gets the dimension of the solution support.
          */
         virtual int getDimension() const
         {
            return getHandle()->dim;
         }

         /**
          * @returns Number of points of the solution.
          */
         virtual int count() const
         {
            return getHandle()->np;
         }

         /**
          * @returns Number of solutions per entity.
          */
         virtual int size() const
         {
            return getHandle()->size;
         }
   };
}

#endif
