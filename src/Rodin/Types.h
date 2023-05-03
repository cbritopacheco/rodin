#ifndef RODIN_TYPES_H
#define RODIN_TYPES_H

#include <cstddef>

#include <boost/unordered_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>

#include <Eigen/Core>

#include "Rodin/Configure.h"

namespace Rodin
{
  using Integer = int;
  using Boolean = bool;
  using Float = double;
  using Scalar = Float;
  using Index = std::size_t;

  template <class T>
  using FlatSet = boost::container::flat_set<T>;

  template <class K, class T>
  using Map = boost::container::map<K, T>;

  template <class K, class T>
  using UnorderedMap = boost::unordered_map<K, T>;

  template <class K, class T>
  using FlatMap = boost::container::flat_map<K, T>;

#if __cpp_size_t_suffix < 202011L
  constexpr
  std::size_t operator "" _UZ (unsigned long long x)
  {
    return x;
  }
#endif
}

#endif
