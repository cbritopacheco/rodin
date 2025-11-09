/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SHARDER_H
#define RODIN_GEOMETRY_SHARDER_H

/**
 * @file
 * @brief Sharder interface for mesh decomposition.
 */

#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Interface for mesh sharding algorithms.
   *
   * A Sharder decomposes a mesh into smaller pieces (shards) that can be
   * processed independently. This is useful for distributed computing and
   * parallel finite element analysis.
   *
   * @tparam Context The context type (e.g., Context::Local, Context::MPI)
   *
   * @see Shard, MeshShard
   */
  template <class Context>
  class Sharder;
}

#endif
