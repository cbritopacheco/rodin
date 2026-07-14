/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Array.h
 * @brief boost::serialization support for Rodin::Array types.
 */
#ifndef RODIN_SERIALIZATION_ARRAY_H
#define RODIN_SERIALIZATION_ARRAY_H

#include <boost/serialization/array.hpp>
#include <boost/serialization/boost_array.hpp>

#include "Rodin/Array.h"

namespace boost::serialization
{
  /**
   * @brief Saves a Rodin array to an archive.
   * @tparam Archive Archive type.
   * @tparam ScalarType Array scalar type.
   * @param ar Archive used for serialization.
   * @param arr Array to save.
   * @param version Serialization version.
   */
  template <class Archive, typename ScalarType>
  void save(Archive& ar, const Rodin::Array<ScalarType>& arr, const unsigned int version)
  {
    (void)version;
    const size_t sz = arr.size();
    ar & sz;
    ar & boost::serialization::make_array(arr.data(), sz);
  }

  /**
   * @brief Loads a Rodin array from an archive.
   * @tparam Archive Archive type.
   * @tparam ScalarType Array scalar type.
   * @param ar Archive used for serialization.
   * @param arr Array to load.
   * @param version Serialization version.
   */
  template <class Archive, typename ScalarType>
  void load(Archive& ar, Rodin::Array<ScalarType>& arr, const unsigned int version)
  {
    (void)version;
    size_t sz;
    ar & sz;
    arr.resize(sz);
    ar & boost::serialization::make_array(arr.data(), sz);
  }

  /**
   * @brief Serializes a Rodin array through split save/load functions.
   * @tparam Archive Archive type.
   * @tparam ScalarType Array scalar type.
   * @param ar Archive used for serialization.
   * @param arr Array to serialize.
   * @param version Serialization version.
   */
  template <class Archive, typename ScalarType>
  void serialize(
      Archive & ar,
      Rodin::Array<ScalarType>& arr,
      const unsigned int version)
  {
    boost::serialization::split_free(ar, arr, version);
  }
}

#endif
