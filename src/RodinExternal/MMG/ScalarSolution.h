/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SCALARSOLUTION_H
#define RODIN_EXTERNAL_MMG_SCALARSOLUTION_H

#include <cassert>

#include <mmg/libmmg.h>
#include <mmg/mmg2d/libmmg2d.h>

#include "Solution.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Scalar solution supported on a mesh.
    *
    * @tparam Dimension Dimension of the solution support
    */
   template <int Dimension>
   class ScalarSolution : public Solution<Dimension, double>
   {};
}

#endif
