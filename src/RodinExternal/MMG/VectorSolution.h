/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_VECTORSOLUTION_H
#define RODIN_EXTERNAL_MMG_VECTORSOLUTION_H

#include <cassert>

#include "Solution.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Vector solution supported on a mesh.
   */
  class VectorSolution : public SolutionBase
  {
    public:
      VectorSolution(MMG5_pSol sol = nullptr)
        : SolutionBase(sol)
      {}

      VectorSolution(const VectorSolution& other)
         : SolutionBase(other)
      {}

      VectorSolution(VectorSolution&& other)
        : SolutionBase(std::move(other))
      {}

      VectorSolution& operator=(VectorSolution&&) = default;
  };
}

#endif

