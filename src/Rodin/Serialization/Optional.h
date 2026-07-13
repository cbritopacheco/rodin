/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Optional.h
 * @brief boost::serialization support for std::optional.
 */
#ifndef RODIN_SERIALIZATION_OPTIONAL_H
#define RODIN_SERIALIZATION_OPTIONAL_H

#include <optional>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost::serialization
{
  /**
   * @brief Saves an optional value to an archive.
   * @tparam Archive Archive type.
   * @tparam T Optional value type.
   * @param ar Archive used for serialization.
   * @param opt Optional value to save.
   * @param version Serialization version.
   */
  template <class Archive, class T>
  void save(Archive& ar, const std::optional<T>& opt, const unsigned int version)
  {
      (void) version;
      bool has_value = opt.has_value();
      ar & has_value;

      if (has_value)
          ar & *opt;
  }

  /**
   * @brief Loads an optional value from an archive.
   * @tparam Archive Archive type.
   * @tparam T Optional value type.
   * @param ar Archive used for serialization.
   * @param opt Optional value to load.
   * @param version Serialization version.
   */
  template <class Archive, class T>
  void load(Archive& ar, std::optional<T>& opt, const unsigned int version)
  {
      (void) version;
      bool has_value;
      ar & has_value;

      if (has_value)
      {
          T value;
          ar & value;
          opt = std::move(value);
      }
      else
      {
          opt.reset();
      }
  }

  /**
   * @brief Serializes an optional value through split save/load functions.
   * @tparam Archive Archive type.
   * @tparam T Optional value type.
   * @param ar Archive used for serialization.
   * @param opt Optional value to serialize.
   * @param version Serialization version.
   */
  template <class Archive, class T>
  void serialize(Archive& ar, std::optional<T>& opt, const unsigned int version)
  {
      split_free(ar, opt, version);
  }
}

#endif
