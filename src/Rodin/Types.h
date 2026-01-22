/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TYPES_H
#define RODIN_TYPES_H

/**
 * @file
 * @brief Fundamental type definitions for the Rodin library.
 *
 * This header defines the standard type aliases used throughout Rodin,
 * including numeric types (Integer, Real, Complex), container types
 * (List, Map, Set), and utility types (Optional, Index).
 *
 * These definitions ensure consistent type usage across the library and
 * facilitate easy customization of underlying types if needed.
 *
 * @see @ref RodinTypes
 */

#include <stack>
#include <cstddef>
#include <complex>
#include <bitset>
#include <optional>

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
  /**
   * @defgroup RodinTypes Type Definitions
   * @brief Standard type aliases and definitions used throughout Rodin
   *
   * This module contains type definitions that provide consistent naming 
   * and standardized types across the entire Rodin library.
   */

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

  /// Standard type for representing lists.
  template <class T>
  using List = boost::container::list<T>;

  /// Standard type for representing deques.
  template <class T>
  using Deque = boost::container::deque<T>;

  /// Standard type for representing stacks.
  /// @ingroup RodinTypes
  template <class T, class Container = Deque<T>>
  using Stack = std::stack<T, Container>;

  /// Standard flat set type.
  /// @ingroup RodinTypes
  template <class T>
  using FlatSet = boost::container::flat_set<T>;

  /// Standard unordered set type.
  /// @ingroup RodinTypes
  template <class T>
  using UnorderedSet = boost::unordered_set<T>;

  /// Standard ordered map type.
  /// @ingroup RodinTypes
  template <class K, class T>
  using Map = boost::container::map<K, T>;

  /// Standard unordered map type.
  /// @ingroup RodinTypes
  template <class ... Params>
  using UnorderedMap = boost::unordered_map<Params...>;

  /// Standard flat map type.
  /// @ingroup RodinTypes
  template <class K, class T>
  using FlatMap = boost::container::flat_map<K, T>;

  /// Standard set of indices.
  /// @ingroup RodinTypes
  using IndexSet = FlatSet<Index>;

  /// Standard vector of indices.
  /// @ingroup RodinTypes
  using IndexVector = std::vector<Index>;

  /// Standard map of indices to arbitrary types.
  /// @ingroup RodinTypes
  template <class T>
  using IndexMap = FlatMap<Index, T>;

  /// Standard bitset template with fixed size.
  /// @ingroup RodinTypes
  template <size_t Size>
  using BitSet = std::bitset<Size>;

  /// Standard bitset with 2 bits.
  /// @ingroup RodinTypes
  using BitSet2 = std::bitset<2>;

  /// Standard optional type wrapper.
  /// @ingroup RodinTypes
  template <class T>
  using Optional = std::optional<T>;

#if __cpp_size_t_suffix < 202011L
  /**
   * @brief User-defined literal for size_t values.
   * @param x Unsigned long long value to convert
   * @return size_t value
   * @note This provides backward compatibility for the _UZ suffix
   */
  constexpr
  std::size_t operator""_UZ (unsigned long long x)
  {
    return x;
  }
#endif
}

#endif
