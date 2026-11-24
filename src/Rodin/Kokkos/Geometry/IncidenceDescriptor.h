/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_GEOMETRY_INCIDENCEDESCRIPTOR_H
#define RODIN_KOKKOS_GEOMETRY_INCIDENCEDESCRIPTOR_H

/**
 * @file
 * @brief Borrowed compressed arrays of one incidence relation.
 */

#include <cstdint>
#include "Rodin/Types.h"

namespace Rodin::Kokkos::Geometry
{

  /**
   * @brief Borrowed resident arrays for one dimension-pair relation.
   * Independent CSR allocations let compute() replace a single relation.
   */
  struct IncidenceDescriptor
  {
    std::uint8_t present; ///< Whether the relation is materialized.
    std::uint8_t prepared; ///< Whether the materialized relation is current.
    const Index* offsets; ///< Row offsets, including the terminal offset.
    const Index* indices; ///< Dimension-local target indices.
    Index count; ///< Number of source rows.
  };

}

#endif
