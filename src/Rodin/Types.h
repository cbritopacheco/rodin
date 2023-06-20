#ifndef RODIN_TYPES_H
#define RODIN_TYPES_H

#include <cstddef>

#include <boost/unordered_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>

#include <Eigen/Core>

#include "Rodin/Configure.h"

#include "Array.h"

namespace Rodin
{
  /// Standard type for representing integer values.
  using Integer = int;

  /// Standard type for representing boolean values.
  using Boolean = bool;

  /// Standard type for representing floating-point values.
  using Float = double;

  /// Standard type for representing scalar values.
  using Scalar = Float;

  /// Standard type for representing indices.
  using Index = std::size_t;

  /// Standard flat set type.
  template <class T>
  using FlatSet = boost::container::flat_set<T>;

  /// Standard ordered map type.
  template <class K, class T>
  using Map = boost::container::map<K, T>;

  /// Standard unordered map type.
  template <class K, class T>
  using UnorderedMap = boost::unordered_map<K, T>;

  /// Standard flat map type.
  template <class K, class T>
  using FlatMap = boost::container::flat_map<K, T>;

  /// Standard set of indices.
  using IndexSet = FlatSet<Index>;

  /// Standard set of indices.
  using IndexArray = Array<Index>;

  /// Standard map of indices.
  template <class T>
  using IndexMap = FlatMap<Index, T>;

#if __cpp_size_t_suffix < 202011L
  constexpr
  std::size_t operator "" _UZ (unsigned long long x)
  {
    return x;
  }
#endif
}

#endif
