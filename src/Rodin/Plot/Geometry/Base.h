/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_H
#define RODIN_CORE_GEOMETRY_H

#include <utility>
#include <type_traits>

#include "Rodin/Math/Common.h"

#include "Concepts.h"

namespace Rodin::Plot::Geometry
{
  template <class T, class Derived>
  class Base
  {
   public:
    // template <class DerivedGeometry>
    // inline
    // constexpr
    // auto intersect(const Geometry<T, DerivedGeometry>& other) const
    //  -> decltype(
    //     std::declval<Derived>().intersect(
    //      std::declval<const Geometry<T, DerivedGeometry>&>()))
    // {
    //  return Derived::intersect(other);
    // }

    // template <class DerivedGeometry>
    // inline
    // constexpr
    // std::optional<LineSegment2D<T>> connect(
    //    const Geometry<T, DerivedGeometry>& other) const
    // {
    //  static_assert(
    //     has_connect_method<Derived, DerivedGeometry>::value,
    //     "Derived::connect(const OtherGeometry&) has not been implemented!"
    //     );
    //  return Derived::connect(other);
    // }

    // template <class DerivedGeometry>
    // inline
    // constexpr
    // typename std::enable_if<
    //  has_distance_method<Derived, DerivedGeometry>::value
    //  && !has_connect_method<Derived, DerivedGeometry>::value,
    // T>::type distance(const Geometry<T, DerivedGeometry>& other) const
    // {
    //  return Derived::distance(other);
    // }

    // template <class DerivedGeometry>
    // inline
    // constexpr
    // typename std::enable_if<
    //  has_connect_method<Derived, DerivedGeometry>::value
    //  && !has_distance_method<Derived, DerivedGeometry>::value,
    // T>::type distance(const Geometry<T, DerivedGeometry>& other) const
    // {
    //  if (constexpr auto lineSegment = Derived::connect(other))
    //    return lineSegment.length();
    //  else
    //    return 0;
    // }

    virtual ~Base() = default;
  };
}

#endif
