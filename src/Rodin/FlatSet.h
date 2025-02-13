/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FLATSET_H
#define RODIN_FLATSET_H

#include "Types.h"

namespace Rodin
{
  /**
   * @brief Functor for comparing two IndexSet objects for equality.
   *
   * This functor defines an operator() that compares two IndexSet objects.
   * Two sets are considered equal if they have the same size and all corresponding
   * elements are equal.
   */
  struct IndexSetEquality
  {
    /**
     * @brief Compares two IndexSet objects for equality.
     *
     * @param lhs The left-hand side IndexSet.
     * @param rhs The right-hand side IndexSet.
     * @return true if both IndexSet objects are equal, false otherwise.
     */
    bool operator()(const IndexSet& lhs, const IndexSet& rhs) const
    {
      if (lhs.size() != rhs.size())
        return false;
      else
        return std::equal(lhs.begin(), lhs.end(), rhs.begin());
    }
  };

  /**
   * @brief Functor for computing a hash value for an IndexSet.
   *
   * This functor defines an operator() that computes a hash value for an IndexSet
   * by combining the hash values of its elements.
   */
  struct IndexSetHash
  {
    /**
     * @brief Computes the hash value for the given IndexSet.
     *
     * @param arr The IndexSet for which to compute the hash.
     * @return A size_t hash value representing the IndexSet.
     */
    size_t operator()(const IndexSet& arr) const
    {
      size_t seed = 0;
      std::for_each(arr.begin(), arr.end(), [&](Rodin::Index v) { boost::hash_combine(seed, v); } );
      return seed;
    }
  };
}

#endif
