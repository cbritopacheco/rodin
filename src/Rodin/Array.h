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
  template <class Real>
  using Array = Eigen::ArrayX<Real>;

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
      assert(lhs.size() > 0);
      assert(rhs.size() > 0);
      if (lhs.size() != rhs.size())
      {
        return false;
      }
      else
      {
        assert(lhs.size() == rhs.size());
        switch (lhs.size())
        {
          case 1:
          {
            return lhs.coeff(0) == rhs.coeff(0);
          }
          case 2:
          {
            return (lhs.coeff(0) == rhs.coeff(0) && lhs.coeff(1) == rhs.coeff(1))
              || (lhs.coeff(0) == rhs.coeff(1) && lhs.coeff(1) == rhs.coeff(0));
          }
          default:
          {
            return std::is_permutation(lhs.begin(), lhs.end(), rhs.begin());
          }
        }
      }
    }
  };

  struct IndexArrayHash
  {
    inline
    size_t operator()(const IndexArray& arr) const
    {
      size_t seed = 0;
      switch (arr.size())
      {
        case 1:
        {
          boost::hash_combine(seed, arr.coeff(0));
          break;
        }
        case 2:
        {
          boost::hash_combine(seed, arr.coeff(0));
          boost::hash_combine(seed, arr.coeff(1));
          break;
        }
        case 3:
        {
          boost::hash_combine(seed, arr.coeff(0));
          boost::hash_combine(seed, arr.coeff(1));
          boost::hash_combine(seed, arr.coeff(2));
          break;
        }
        case 4:
        {
          boost::hash_combine(seed, arr.coeff(0));
          boost::hash_combine(seed, arr.coeff(1));
          boost::hash_combine(seed, arr.coeff(2));
          boost::hash_combine(seed, arr.coeff(3));
          break;
        }
        case 5:
        {
          boost::hash_combine(seed, arr.coeff(0));
          boost::hash_combine(seed, arr.coeff(1));
          boost::hash_combine(seed, arr.coeff(2));
          boost::hash_combine(seed, arr.coeff(3));
          boost::hash_combine(seed, arr.coeff(4));
          break;
        }
        default:
        {
          std::for_each(arr.begin(), arr.end(), [&](Rodin::Index v) { boost::hash_combine(seed, v); } );
          break;
        }
      }
      return seed;
    }
  };

  struct IndexArraySymmetricHash
  {
    inline
    size_t operator()(const IndexArray& arr) const
    {
      size_t seed = 0;
      switch (arr.size())
      {
        case 0:
        {
          break;
        }
        case 1:
        {
          seed = boost::hash_value(arr.coeff(0));
          break;
        }
        case 2:
        {
          seed = boost::hash_value(arr.coeff(0))
               + boost::hash_value(arr.coeff(1));
          break;
        }
        case 3:
        {
          seed = boost::hash_value(arr.coeff(0))
               + boost::hash_value(arr.coeff(1))
               + boost::hash_value(arr.coeff(2));
          break;
        }
        case 4:
        {
          seed = boost::hash_value(arr.coeff(0))
               + boost::hash_value(arr.coeff(1))
               + boost::hash_value(arr.coeff(2))
               + boost::hash_value(arr.coeff(3));
          break;
        }
        default:
        {
          std::for_each(arr.begin(), arr.end(), [&](auto&& v) { seed += boost::hash_value(v); } );
          break;
        }
      }
      return seed;
    }
  };
}

#endif

