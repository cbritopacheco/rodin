/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_CUDSSVECTOR_H
#define RODIN_KOKKOS_SOLVER_CUDSSVECTOR_H

/**
 * @file
 * @brief Direct-solver access to backend vector storage.
 */

#include "Rodin/Solver/CuDSSOperator.h"
#include "Rodin/Kokkos/Math/Vector.h"

namespace Rodin::Solver
{
  /// @brief Backend vectors are read directly when host access is unavailable.
  template <class Scalar, class Backend>
  struct CuDSSVector<::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
  {
    /// @brief Whether the allocation requires backend access.
    static constexpr bool Resident = !::Kokkos::SpaceAccessibility<::Kokkos::HostSpace,
      typename ::Rodin::Kokkos::Math::Vector<Scalar, Backend>::MemorySpace>::accessible;
  };
}

#endif
