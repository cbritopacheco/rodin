/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_CONCEPTS_H
#define RODIN_CORE_GEOMETRY_CONCEPTS_H

#include <utility>
#include <type_traits>

namespace Rodin::Plot::Geometry
{
  template<class T, class OtherGeometry, class = void>
  struct has_intersect_method : std::false_type{};

  template<class T, class OtherGeometry>
  struct has_intersect_method<T, OtherGeometry,
   std::void_t<decltype(std::declval<T>().intersect(
       std::declval<const OtherGeometry&>()))>>
    : std::true_type{};

  template<class T, class OtherGeometry, class = void>
  struct has_distance_method : std::false_type{};

  template<class T, class OtherGeometry>
  struct has_distance_method<T, OtherGeometry,
   std::void_t<decltype(std::declval<T>().distance(
       std::declval<const OtherGeometry&>()))>>
    : std::true_type{};

  template<class T, class OtherGeometry, class = void>
  struct has_connect_method : std::false_type{};

  template<class T, class OtherGeometry>
  struct has_connect_method<T, OtherGeometry,
   std::void_t<decltype(std::declval<T>().connect(
       std::declval<const OtherGeometry&>()))>>
    : std::true_type{};
}
#endif
