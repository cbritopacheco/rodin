/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_CONTEXT_HIP_H
#define RODIN_KOKKOS_CONTEXT_HIP_H

#include "Rodin/Context/ForwardDecls.h"

namespace Rodin::Context
{
  /**
   * @brief Selects HIP execution in a composed Kokkos context.
   *
   * The tag carries no allocation, runtime initialization, or distribution
   * state. @ref Kokkos maps it to the HIP execution space supplied by Kokkos.
   */
  class HIP final
  {};
}

#endif
