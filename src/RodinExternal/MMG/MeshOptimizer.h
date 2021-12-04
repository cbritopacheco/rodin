/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHOPTIMIZER_H
#define RODIN_MESH_MESHOPTIMIZER_H

#include "ScalarSolution2D.h"

namespace Rodin::External::MMG
{
  template <class MeshType, class Derived>
  class MeshOptimizer
  {
    public:
      ScalarSolution2D optimize(MeshType& mesh)
      {
        return Derived::optimize(mesh);
      }
  };
}

#endif
