/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PolytopeType.h
 * @brief boost::serialization support for Geometry::Polytope::Type.
 */
#ifndef RODIN_SERIALIZATION_POLYTOPETYPE_H
#define RODIN_SERIALIZATION_POLYTOPETYPE_H

#include <type_traits>
#include <boost/serialization/array.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "Rodin/Array.h"

namespace boost::serialization
{
  /**
   * @brief Serializes a polytope type enum value.
   * @tparam Archive Archive type.
   * @tparam ScalarType Unused scalar tag kept for compatibility.
   * @param ar Archive used for serialization.
   * @param t Polytope type to serialize.
   * @param version Serialization version.
   */
  template <class Archive, class ScalarType>
  void serialize(
      Archive & ar, const Rodin::Geometry::Polytope::Type& t, const unsigned int version)
  {
    (void) version;
    ar & static_cast<std::underlying_type_t<Rodin::Geometry::Polytope::Type>>(t);
  }
}

#endif

