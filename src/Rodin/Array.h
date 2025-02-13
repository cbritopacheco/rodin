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

  /**
   * @brief Functor for comparing two index arrays for symmetric equality.
   *
   * This functor compares two index arrays and considers them equal if they
   * contain the same elements, regardless of their order. For arrays of size 1
   * and 2, special cases are handled explicitly for efficiency.
   * For larger arrays, std::is_permutation is used.
   */
  struct IndexArraySymmetricEquality
  {
    /**
     * @brief Compares two index arrays for symmetric equality.
     *
     * @param lhs The left-hand side index array.
     * @param rhs The right-hand side index array.
     * @return true if the arrays are symmetric equal, false otherwise.
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

  /**
   * @brief Functor for computing a hash value for an index array.
   *
   * This functor computes a hash value for an index array using boost::hash_combine.
   * Special cases are handled for arrays of sizes 1 to 5 for efficiency, while larger arrays
   * are processed using a generic loop.
   */
  struct IndexArrayHash
  {
    /**
     * @brief Computes the hash value for the given index array.
     *
     * @param arr The index array for which to compute the hash.
     * @return A size_t hash value representing the index array.
     */
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

  /**
   * @brief Functor for computing a symmetric hash value for an index array.
   *
   * This functor computes a hash value for an index array in a way that is
   * independent of the order of the elements. It does so by summing the hash
   * values of individual elements using boost::hash_value.
   */
  struct IndexArraySymmetricHash
  {
    /**
     * @brief Computes the symmetric hash value for the given index array.
     *
     * @param arr The index array for which to compute the symmetric hash.
     * @return A size_t hash value representing the index array independent of element order.
     */
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

