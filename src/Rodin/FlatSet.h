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
  struct IndexSetEquality
  {
    inline
    bool operator()(const IndexSet& lhs, const IndexSet& rhs) const
    {
      if (lhs.size() != rhs.size())
        return false;
      else
        return std::equal(lhs.begin(), lhs.end(), rhs.begin());
    }
  };

  struct IndexSetHash
  {
    inline
    size_t operator()(const IndexSet& arr) const
    {
      size_t seed = 0;
      std::for_each(arr.begin(), arr.end(), [&](Rodin::Index v) { boost::hash_combine(seed, v); } );
      return seed;
    }
  };
}
#endif
