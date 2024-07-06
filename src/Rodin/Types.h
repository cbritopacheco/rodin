#ifndef RODIN_TYPES_H
#define RODIN_TYPES_H

#include <stack>
#include <cstddef>
#include <complex>

#include <boost/unordered_map.hpp>
#include <boost/container/map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/deque.hpp>
#include <boost/container/list.hpp>

#include <Eigen/Core>

#include "Rodin/Configure.h"

namespace Rodin
{
  /// Standard type for representing integer values.
  using Integer = int;

  /// Standard type for representing boolean values.
  using Boolean = bool;

  /// Standard type for representing (32 bit) floating-point values.
  using Float = float;

  /// Standard type for representing double precision (64 bit) floating-point values.
  using Double = double;

  /// Standard type for representing indices.
  using Index = std::size_t;

  /// Standard type for representing scalar values.
  using Real = Double;

  /// Standard type for representing complex values.
  using Complex = std::complex<Real>;

  template <class T>
  using List = boost::container::list<T>;

  template <class T>
  using Deque = boost::container::deque<T>;

  template <class T, class Container = Deque<T>>
  using Stack = std::stack<T, Container>;

  /// Standard flat set type.
  template <class T>
  using FlatSet = boost::container::flat_set<T>;

  /// Standard unordered set type.
  template <class T>
  using UnorderedSet = boost::unordered_set<T>;

  /// Standard ordered map type.
  template <class K, class T>
  using Map = boost::container::map<K, T>;

  /// Standard unordered map type.
  template <class ... Params>
  using UnorderedMap = boost::unordered_map<Params...>;

  /// Standard flat map type.
  template <class K, class T>
  using FlatMap = boost::container::flat_map<K, T>;

  /// Standard set of indices.
  using IndexSet = FlatSet<Index>;

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
