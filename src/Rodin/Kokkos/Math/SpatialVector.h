/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MATH_SPATIALVECTOR_H
#define RODIN_KOKKOS_MATH_SPATIALVECTOR_H

/**
 * @file
 * @brief Device-callable bounded spatial vector.
 */

#include <Kokkos_Core.hpp>
#include "Rodin/Types.h"

namespace Rodin::Kokkos::Math
{
  /**
   * @brief Device-local vector with at most three active coordinates.
   *
   * Storage is independent of host Eigen evaluation and never allocates.
   */
  class SpatialVector
  {
    public:
      /// @brief Constructs the empty coordinate tuple of a reference point.
      KOKKOS_INLINE_FUNCTION SpatialVector()
        : m_size(0),
          m_data{}
      {}

      /**
       * @brief Constructs a zero vector.
       * @param size Active dimension, between zero and three.
       */
      KOKKOS_INLINE_FUNCTION explicit SpatialVector(Index size)
        : m_size(size),
          m_data{}
      {
        KOKKOS_ASSERT(size <= 3);
      }

      /// @brief Returns the active dimension. @return Coordinate count.
      KOKKOS_INLINE_FUNCTION Index size() const
      {
        return m_size;
      }

      /// @brief Returns the active dimension. @return Coordinate count.
      KOKKOS_INLINE_FUNCTION Index getDimension() const
      {
        return m_size;
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Mutable coordinate.
       */
      KOKKOS_INLINE_FUNCTION Real& operator()(Index i)
      {
        KOKKOS_ASSERT(i < m_size);
        return m_data[i];
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Coordinate value.
       */
      KOKKOS_INLINE_FUNCTION Real operator()(Index i) const
      {
        KOKKOS_ASSERT(i < m_size);
        return m_data[i];
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Mutable coordinate.
       */
      KOKKOS_INLINE_FUNCTION Real& operator[](Index i)
      {
        return operator()(i);
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Coordinate value.
       */
      KOKKOS_INLINE_FUNCTION Real operator[](Index i) const
      {
        return operator()(i);
      }

    private:
      Index m_size;
      ::Kokkos::Array<Real, 3> m_data;
  };
}

#endif
