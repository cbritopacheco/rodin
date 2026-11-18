/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_INTEGRATION_H
#define RODIN_KOKKOS_VARIATIONAL_INTEGRATION_H

/**
 * @file
 * @brief Compatibility umbrella for contribution containers and their layout.
 */

#include "Rodin/Kokkos/Math/Triplets.h"
#include "Rodin/Kokkos/Math/Entries.h"
#include "Layout.h"

namespace Rodin::Kokkos::Variational
{
  /// @brief Compatibility name of sparse matrix contributions.
  template <class MemorySpace>
  using Triplets = ::Rodin::Kokkos::Math::Triplets<MemorySpace>;

  /// @brief Compatibility name of sparse vector contributions.
  template <class MemorySpace>
  using Entries = ::Rodin::Kokkos::Math::Entries<MemorySpace>;
}

#endif
