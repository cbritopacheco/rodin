/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_CONCEPTS_H
#define RODIN_GEOMETRY_EUCLIDEAN_CONCEPTS_H

/**
 * @file
 * @brief Type traits for detecting geometric operations on Euclidean objects.
 */

#include <utility>
#include <type_traits>

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Type trait to detect if a type has an intersect method.
   *
   * Detects whether type @p T has a member function `intersect(const OtherGeometry&)`.
   *
   * @tparam T Type to check
   * @tparam OtherGeometry Argument type for intersect method
   *
   * Provides the member constant `value` which is equal to `true` if @p T has
   * an `intersect` method accepting @p OtherGeometry, or `false` otherwise.
   */
  template<class T, class OtherGeometry, class = void>
  struct has_intersect_method : std::false_type{};

  /**
   * @brief Specialization when intersect method exists.
   */
  template<class T, class OtherGeometry>
  struct has_intersect_method<T, OtherGeometry,
   std::void_t<decltype(std::declval<T>().intersect(
       std::declval<const OtherGeometry&>()))>>
    : std::true_type{};

  /**
   * @brief Type trait to detect if a type has a distance method.
   *
   * Detects whether type @p T has a member function `distance(const OtherGeometry&)`.
   *
   * @tparam T Type to check
   * @tparam OtherGeometry Argument type for distance method
   *
   * Provides the member constant `value` which is equal to `true` if @p T has
   * a `distance` method accepting @p OtherGeometry, or `false` otherwise.
   */
  template<class T, class OtherGeometry, class = void>
  struct has_distance_method : std::false_type{};

  /**
   * @brief Specialization when distance method exists.
   */
  template<class T, class OtherGeometry>
  struct has_distance_method<T, OtherGeometry,
   std::void_t<decltype(std::declval<T>().distance(
       std::declval<const OtherGeometry&>()))>>
    : std::true_type{};

  /**
   * @brief Type trait to detect if a type has a connect method.
   *
   * Detects whether type @p T has a member function `connect(const OtherGeometry&)`.
   *
   * @tparam T Type to check
   * @tparam OtherGeometry Argument type for connect method
   *
   * Provides the member constant `value` which is equal to `true` if @p T has
   * a `connect` method accepting @p OtherGeometry, or `false` otherwise.
   */
  template<class T, class OtherGeometry, class = void>
  struct has_connect_method : std::false_type{};

  /**
   * @brief Specialization when connect method exists.
   */
  template<class T, class OtherGeometry>
  struct has_connect_method<T, OtherGeometry,
   std::void_t<decltype(std::declval<T>().connect(
       std::declval<const OtherGeometry&>()))>>
    : std::true_type{};
}
#endif
