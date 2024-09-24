/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CONTEXT_SEQUENTIAL_H
#define RODIN_CONTEXT_SEQUENTIAL_H

#include "Base.h"

namespace Rodin::Context
{
  /**
   * @brief Represents a single machine context.
   *
   * The Local context refers to an execution model where operations are
   * confined to a single machine or node, utilizing shared memory without the
   * need for distributed computing. While operating within a local scope, this
   * context can leverage multithreading or parallelism, but it does not
   * involve communication across multiple machines.
   */
  class Local : public Base
  {};
}

#endif


