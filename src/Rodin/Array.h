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
    static inline Rodin::Index get(const Rodin::IndexArray& a, size_t i)
    {
      return a.coeff(static_cast<Eigen::Index>(i));
    }

    static inline void cswap(Rodin::Index& x, Rodin::Index& y)
    {
      if (y < x) std::swap(x, y);
    }

    // Branch-light sort networks (ascending)
    static inline void sort3(std::array<Rodin::Index,3>& a)
    {
      cswap(a[0], a[1]);
      cswap(a[1], a[2]);
      cswap(a[0], a[1]);
    }

    static inline void sort4(std::array<Rodin::Index,4>& a)
    {
      cswap(a[0], a[1]);
      cswap(a[2], a[3]);
      cswap(a[0], a[2]);
      cswap(a[1], a[3]);
      cswap(a[1], a[2]);
    }

    // 6-input sorting network (12 compare-exchanges)
    static inline void sort6(std::array<Rodin::Index,6>& a)
    {
      cswap(a[1], a[2]);
      cswap(a[4], a[5]);
      cswap(a[0], a[2]);
      cswap(a[3], a[5]);
      cswap(a[0], a[1]);
      cswap(a[3], a[4]);
      cswap(a[2], a[5]);
      cswap(a[0], a[3]);
      cswap(a[1], a[4]);
      cswap(a[2], a[4]);
      cswap(a[1], a[3]);
      cswap(a[2], a[3]);
    }

    // 8-input sorting network (19 compare-exchanges; good practical choice)
    static inline void sort8(std::array<Rodin::Index,8>& a)
    {
      cswap(a[0], a[1]); cswap(a[2], a[3]); cswap(a[4], a[5]); cswap(a[6], a[7]);
      cswap(a[0], a[2]); cswap(a[1], a[3]); cswap(a[4], a[6]); cswap(a[5], a[7]);
      cswap(a[1], a[2]); cswap(a[5], a[6]);
      cswap(a[0], a[4]); cswap(a[1], a[5]); cswap(a[2], a[6]); cswap(a[3], a[7]);
      cswap(a[2], a[4]); cswap(a[3], a[5]);
      cswap(a[1], a[2]); cswap(a[3], a[4]); cswap(a[5], a[6]);
    }

    // Exact unrolled multiset-equality for size 3 (6 cases)
    static inline bool eq3(const Rodin::IndexArray& lhs, const Rodin::IndexArray& rhs)
    {
      const auto a0 = get(lhs,0), a1 = get(lhs,1), a2 = get(lhs,2);
      const auto b0 = get(rhs,0), b1 = get(rhs,1), b2 = get(rhs,2);

      return (a0==b0 && a1==b1 && a2==b2) ||
             (a0==b0 && a1==b2 && a2==b1) ||
             (a0==b1 && a1==b0 && a2==b2) ||
             (a0==b1 && a1==b2 && a2==b0) ||
             (a0==b2 && a1==b0 && a2==b1) ||
             (a0==b2 && a1==b1 && a2==b0);
    }

    // Exact unrolled multiset-equality for size 4 using sorting network (faster than 24 perms)
    static inline bool eq4(const Rodin::IndexArray& lhs, const Rodin::IndexArray& rhs)
    {
      std::array<Rodin::Index,4> a{ get(lhs,0), get(lhs,1), get(lhs,2), get(lhs,3) };
      std::array<Rodin::Index,4> b{ get(rhs,0), get(rhs,1), get(rhs,2), get(rhs,3) };
      sort4(a);
      sort4(b);
      return a == b;
    }

    bool operator()(const Rodin::IndexArray& lhs, const Rodin::IndexArray& rhs) const
    {
      if (lhs.size() != rhs.size())
        return false;

      const size_t n = static_cast<size_t>(lhs.size());

      switch (n)
      {
        case 0: return true;
        case 1: return lhs.coeff(0) == rhs.coeff(0);
        case 2:
          return (lhs.coeff(0) == rhs.coeff(0) && lhs.coeff(1) == rhs.coeff(1))
              || (lhs.coeff(0) == rhs.coeff(1) && lhs.coeff(1) == rhs.coeff(0));
        case 3:
          // Often fastest for triangles: no temp arrays, just comparisons
          return eq3(lhs, rhs);
        case 4:
          return eq4(lhs, rhs);
        case 6:
        {
          std::array<Rodin::Index,6> a{ get(lhs,0), get(lhs,1), get(lhs,2), get(lhs,3), get(lhs,4), get(lhs,5) };
          std::array<Rodin::Index,6> b{ get(rhs,0), get(rhs,1), get(rhs,2), get(rhs,3), get(rhs,4), get(rhs,5) };
          sort6(a);
          sort6(b);
          return a == b;
        }
        case 8:
        {
          std::array<Rodin::Index,8> a{ get(lhs,0), get(lhs,1), get(lhs,2), get(lhs,3), get(lhs,4), get(lhs,5), get(lhs,6), get(lhs,7) };
          std::array<Rodin::Index,8> b{ get(rhs,0), get(rhs,1), get(rhs,2), get(rhs,3), get(rhs,4), get(rhs,5), get(rhs,6), get(rhs,7) };
          sort8(a);
          sort8(b);
          return a == b;
        }
        default:
        {
          // Rare path
          std::vector<Rodin::Index> a(n), b(n);
          for (size_t i = 0; i < n; ++i)
          {
            a[i] = get(lhs, i);
            b[i] = get(rhs, i);
          }
          std::sort(a.begin(), a.end());
          std::sort(b.begin(), b.end());
          return a == b;
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
    static inline uint64_t splitmix64(uint64_t x)
    {
      x += 0x9e3779b97f4a7c15ULL;
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
      x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
      return x ^ (x >> 31);
    }

    static inline uint64_t rotl64(uint64_t x, int r)
    {
      return (x << r) | (x >> (64 - r));
    }

    static inline size_t fold_to_size_t(uint64_t x)
    {
      if constexpr (sizeof(size_t) >= sizeof(uint64_t))
        return static_cast<size_t>(x);
      else
        return static_cast<size_t>(x ^ (x >> 32)); // fold to 32-bit
    }

    /**
     * @brief Computes the symmetric hash value for the given index array.
     *
     * @param arr The index array for which to compute the symmetric hash.
     * @return A size_t hash value representing the index array independent of element order.
     */
    size_t operator()(const Rodin::IndexArray& arr) const
    {
      const auto n = static_cast<size_t>(arr.size());

      // Fast paths for common entity sizes
      switch (n)
      {
        case 0:
        {
          return fold_to_size_t(splitmix64(0x243f6a8885a308d3ULL));
        }
        case 1:
        {
          const uint64_t h0 = splitmix64(static_cast<uint64_t>(arr.coeff(0)));
          uint64_t x = h0 ^ splitmix64(1);
          return fold_to_size_t(splitmix64(x));
        }
        case 2:
        {
          const uint64_t h0 = splitmix64(static_cast<uint64_t>(arr.coeff(0)));
          const uint64_t h1 = splitmix64(static_cast<uint64_t>(arr.coeff(1)));

          uint64_t a = h0 + h1;
          uint64_t b = rotl64(h0, 23) ^ rotl64(h1, 23);
          uint64_t c = (h0 * (h0 | 1ULL)) + (h1 * (h1 | 1ULL));
          uint64_t d = (h0 | 1ULL) * (h1 | 1ULL);

          uint64_t x = a ^ rotl64(b, 17) ^ rotl64(c, 31) ^ rotl64(d, 47) ^ splitmix64(2);
          return fold_to_size_t(splitmix64(x));
        }
        case 3:
        {
          const uint64_t h0 = splitmix64(static_cast<uint64_t>(arr.coeff(0)));
          const uint64_t h1 = splitmix64(static_cast<uint64_t>(arr.coeff(1)));
          const uint64_t h2 = splitmix64(static_cast<uint64_t>(arr.coeff(2)));

          uint64_t a = h0 + h1 + h2;
          uint64_t b = rotl64(h0, 23) ^ rotl64(h1, 23) ^ rotl64(h2, 23);
          uint64_t c = (h0 * (h0 | 1ULL)) + (h1 * (h1 | 1ULL)) + (h2 * (h2 | 1ULL));
          uint64_t d = (h0 | 1ULL) * (h1 | 1ULL) * (h2 | 1ULL);

          uint64_t x = a ^ rotl64(b, 17) ^ rotl64(c, 31) ^ rotl64(d, 47) ^ splitmix64(3);
          return fold_to_size_t(splitmix64(x));
        }
        case 4:
        {
          const uint64_t h0 = splitmix64(static_cast<uint64_t>(arr.coeff(0)));
          const uint64_t h1 = splitmix64(static_cast<uint64_t>(arr.coeff(1)));
          const uint64_t h2 = splitmix64(static_cast<uint64_t>(arr.coeff(2)));
          const uint64_t h3 = splitmix64(static_cast<uint64_t>(arr.coeff(3)));

          uint64_t a = h0 + h1 + h2 + h3;
          uint64_t b = rotl64(h0, 23) ^ rotl64(h1, 23) ^ rotl64(h2, 23) ^ rotl64(h3, 23);
          uint64_t c = (h0 * (h0 | 1ULL)) + (h1 * (h1 | 1ULL)) + (h2 * (h2 | 1ULL)) + (h3 * (h3 | 1ULL));
          uint64_t d = (h0 | 1ULL) * (h1 | 1ULL) * (h2 | 1ULL) * (h3 | 1ULL);

          uint64_t x = a ^ rotl64(b, 17) ^ rotl64(c, 31) ^ rotl64(d, 47) ^ splitmix64(4);
          return fold_to_size_t(splitmix64(x));
        }
        case 6:
        {
          const uint64_t h0 = splitmix64(static_cast<uint64_t>(arr.coeff(0)));
          const uint64_t h1 = splitmix64(static_cast<uint64_t>(arr.coeff(1)));
          const uint64_t h2 = splitmix64(static_cast<uint64_t>(arr.coeff(2)));
          const uint64_t h3 = splitmix64(static_cast<uint64_t>(arr.coeff(3)));
          const uint64_t h4 = splitmix64(static_cast<uint64_t>(arr.coeff(4)));
          const uint64_t h5 = splitmix64(static_cast<uint64_t>(arr.coeff(5)));

          const uint64_t a = h0 + h1 + h2 + h3 + h4 + h5;
          const uint64_t b = rotl64(h0, 23) ^ rotl64(h1, 23) ^ rotl64(h2, 23)
                           ^ rotl64(h3, 23) ^ rotl64(h4, 23) ^ rotl64(h5, 23);
          const uint64_t c = (h0 * (h0 | 1ULL)) + (h1 * (h1 | 1ULL)) + (h2 * (h2 | 1ULL))
                           + (h3 * (h3 | 1ULL)) + (h4 * (h4 | 1ULL)) + (h5 * (h5 | 1ULL));
          const uint64_t d = (h0 | 1ULL) * (h1 | 1ULL) * (h2 | 1ULL)
                           * (h3 | 1ULL) * (h4 | 1ULL) * (h5 | 1ULL);

          const uint64_t x = a ^ rotl64(b, 17) ^ rotl64(c, 31) ^ rotl64(d, 47) ^ splitmix64(6);
          return fold_to_size_t(splitmix64(x));
        }
        case 8:
        {
          const uint64_t h0 = splitmix64(static_cast<uint64_t>(arr.coeff(0)));
          const uint64_t h1 = splitmix64(static_cast<uint64_t>(arr.coeff(1)));
          const uint64_t h2 = splitmix64(static_cast<uint64_t>(arr.coeff(2)));
          const uint64_t h3 = splitmix64(static_cast<uint64_t>(arr.coeff(3)));
          const uint64_t h4 = splitmix64(static_cast<uint64_t>(arr.coeff(4)));
          const uint64_t h5 = splitmix64(static_cast<uint64_t>(arr.coeff(5)));
          const uint64_t h6 = splitmix64(static_cast<uint64_t>(arr.coeff(6)));
          const uint64_t h7 = splitmix64(static_cast<uint64_t>(arr.coeff(7)));

          const uint64_t a = h0 + h1 + h2 + h3 + h4 + h5 + h6 + h7;
          const uint64_t b = rotl64(h0, 23) ^ rotl64(h1, 23) ^ rotl64(h2, 23) ^ rotl64(h3, 23)
                           ^ rotl64(h4, 23) ^ rotl64(h5, 23) ^ rotl64(h6, 23) ^ rotl64(h7, 23);
          const uint64_t c = (h0 * (h0 | 1ULL)) + (h1 * (h1 | 1ULL)) + (h2 * (h2 | 1ULL)) + (h3 * (h3 | 1ULL))
                           + (h4 * (h4 | 1ULL)) + (h5 * (h5 | 1ULL)) + (h6 * (h6 | 1ULL)) + (h7 * (h7 | 1ULL));
          const uint64_t d = (h0 | 1ULL) * (h1 | 1ULL) * (h2 | 1ULL) * (h3 | 1ULL)
                           * (h4 | 1ULL) * (h5 | 1ULL) * (h6 | 1ULL) * (h7 | 1ULL);

          const uint64_t x = a ^ rotl64(b, 17) ^ rotl64(c, 31) ^ rotl64(d, 47) ^ splitmix64(8);
          return fold_to_size_t(splitmix64(x));
        }
        default:
        {
          uint64_t a = 0;
          uint64_t b = 0;
          uint64_t c = 0x9e3779b97f4a7c15ULL;
          uint64_t d = 0xbf58476d1ce4e5b9ULL;

          for (Rodin::Index v : arr)
          {
            const uint64_t h = splitmix64(static_cast<uint64_t>(v));
            a += h;
            b ^= rotl64(h, 23);
            c += h * (h | 1ULL);
            d *= (h | 1ULL);
          }

          uint64_t x = a ^ rotl64(b, 17) ^ rotl64(c, 31) ^ rotl64(d, 47) ^ splitmix64(static_cast<uint64_t>(n));
          return fold_to_size_t(splitmix64(x));
        }
      }
    }
  };
}

#endif

