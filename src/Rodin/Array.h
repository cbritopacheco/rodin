/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ARRAY_H
#define RODIN_ARRAY_H

/**
 * @file
 * @brief Defines array types and utility functors for index arrays.
 *
 * This header provides Eigen-based array type aliases and various functors
 * for comparing, hashing, and manipulating index arrays used throughout
 * the Rodin library.
 */

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Types.h"

namespace Rodin
{
  /**
   * @brief Alias for a dynamically sized array.
   *
   * This template alias defines a standard array type based on Eigen::ArrayX.
   *
   * @tparam ScalarType The scalar type of the array.
   */
  template <class ScalarType>
  using Array = Eigen::ArrayX<ScalarType>;

  /**
   * @brief Alias for an index array.
   *
   * This alias defines an index array using the standard array type with
   * a predefined index type.
   */
  using IndexArray = Array<Index>;

  struct IndexArrayCompare
  {
    bool operator()(const IndexArray& lhs, const IndexArray& rhs) const
    {
      return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
    }
  };

  /**
   * @brief Functor for comparing two index arrays for equality.
   *
   * This functor provides an operator() that returns true if the two given
   * index arrays are equal. Two arrays are considered equal if they are both
   * empty, or if they have the same size and all corresponding elements are
   * equal.
   */
  struct IndexArrayEquality
  {
    /**
     * @brief Compares two index arrays for equality.
     *
     * @param lhs The left-hand side index array.
     * @param rhs The right-hand side index array.
     * @return true if the arrays are equal, false otherwise.
     */
    bool operator()(const IndexArray& lhs, const IndexArray& rhs) const
    {
      if (lhs.size() == 0 && rhs.size() == 0)
      {
        return true;
      }
      else if (lhs.size() != rhs.size())
      {
        return false;
      }
      else
      {
        return (lhs == rhs).all();
      }
    }
  };
}

#endif

