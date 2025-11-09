/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEXCOUNT_H
#define RODIN_GEOMETRY_SIMPLEXCOUNT_H

/**
 * @file
 * @brief Utility class for counting polytopes by dimension.
 */

#include <vector>

namespace Rodin::Geometry
{
  /**
   * @brief Container for counting polytopes in each dimension.
   *
   * This class maintains counts of polytopes for each dimension from 0
   * (vertices) up to the mesh dimension. It provides convenient indexing
   * and initialization for tracking mesh topology.
   */
  class SimplexCount
  {
    public:
      /**
       * @brief Default constructor.
       */
      SimplexCount() = default;

      /**
       * @brief Constructs a SimplexCount for a mesh of given dimension.
       * @param[in] meshDimension The topological dimension of the mesh
       *
       * Initializes counts for dimensions 0 through @p meshDimension to zero.
       */
      SimplexCount(size_t meshDimension)
        : m_counts(meshDimension + 1, 0)
      {}

      /**
       * @brief Constructs from an initializer list.
       * @param[in] l List of counts for each dimension
       */
      SimplexCount(std::initializer_list<size_t> l)
        : m_counts(l)
      {}

      /**
       * @brief Copy constructor.
       */
      SimplexCount(const SimplexCount&) = default;

      /**
       * @brief Move constructor.
       */
      SimplexCount(SimplexCount&&) = default;

      /**
       * @brief Copy assignment operator.
       */
      SimplexCount& operator=(const SimplexCount&) = default;

      /**
       * @brief Move assignment operator.
       */
      SimplexCount& operator=(SimplexCount&&) = default;

      /**
       * @brief Initializes the count structure for a mesh of given dimension.
       * @param[in] meshDimension The topological dimension of the mesh
       * @returns Reference to this object for method chaining
       *
       * Resizes internal storage to hold counts for dimensions 0 through
       * @p meshDimension, initializing all counts to zero.
       */
      inline
      SimplexCount& initialize(size_t meshDimension)
      {
        m_counts.resize(meshDimension + 1, 0);
        return *this;
      }

      /**
       * @brief Accesses the count at dimension @p d (mutable).
       * @param[in] d Dimension index
       * @returns Reference to the count at dimension @p d
       */
      inline
      size_t& at(size_t d)
      {
        return m_counts[d];
      }

      /**
       * @brief Accesses the count at dimension @p d (const).
       * @param[in] d Dimension index
       * @returns Const reference to the count at dimension @p d
       */
      inline
      const size_t& at(size_t d) const
      {
        return m_counts[d];
      }

    private:
      std::vector<size_t> m_counts;
  };
}

#endif
