/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ARRAY_H
#define RODIN_ARRAY_H

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Types.h"

namespace Rodin
{
  /// Standard array type
  template <class Scalar>
  using Array = Eigen::ArrayX<Scalar>;

  /// Standard array of indices.
  using IndexArray = Array<Index>;

  struct IndexArrayEquality
  {
    inline
    bool operator()(const IndexArray& lhs, const IndexArray& rhs) const
    {
      return (lhs == rhs).all();
    }
  };

  struct IndexArraySymmetricEquality
  {
    inline
    bool operator()(const IndexArray& lhs, const IndexArray& rhs) const
    {
      return std::is_permutation(lhs.begin(), lhs.end(), rhs.begin());
    }
  };

  struct IndexArrayHash
  {
    inline
    size_t operator()(const IndexArray& arr) const
    {
      size_t seed = 0;
      std::for_each(arr.begin(), arr.end(), [&](Rodin::Index v) { boost::hash_combine(seed, v); } );
      return seed;
    }
  };

  struct IndexArraySymmetricHash
  {
    inline
    size_t operator()(const IndexArray& arr) const
    {
      size_t seed = 0;
      std::for_each(arr.begin(), arr.end(), [&](auto&& v) { seed += boost::hash_value(v); } );
      return seed;
    }
  };
}

#endif

