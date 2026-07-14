/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FlatMap.h
 * @brief boost::serialization support for boost::container::flat_map.
 */
#ifndef RODIN_SERIALIZATION_FLATMAP_H
#define RODIN_SERIALIZATION_FLATMAP_H

#include <boost/container/flat_map.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/utility.hpp>

namespace boost::serialization
{
  /**
   * @brief Saves a flat map to an archive.
   * @tparam Archive Archive type.
   * @tparam Key Map key type.
   * @tparam T Map mapped value type.
   * @tparam Compare Key comparison type.
   * @tparam Allocator Allocator type.
   * @param ar Archive used for serialization.
   * @param map Flat map to save.
   * @param version Serialization version.
   */
  template <class Archive, typename Key, typename T, typename Compare, typename Allocator>
  void save(Archive& ar,
    const boost::container::flat_map<Key, T, Compare, Allocator>& map,
    const unsigned int version)
  {
    (void)version;
    std::size_t size = map.size();
    ar & size;
    for (const auto& pair : map)
    {
      ar & pair;
    }
  }

  /**
   * @brief Loads a flat map from an archive.
   * @tparam Archive Archive type.
   * @tparam Key Map key type.
   * @tparam T Map mapped value type.
   * @tparam Compare Key comparison type.
   * @tparam Allocator Allocator type.
   * @param ar Archive used for serialization.
   * @param map Flat map to load.
   * @param version Serialization version.
   */
  template <class Archive, typename Key, typename T, typename Compare, typename Allocator>
  void load(Archive& ar, boost::container::flat_map<Key, T, Compare, Allocator>& map,
    const unsigned int version)
  {
    (void)version;
    std::size_t size;
    ar & size;
    map.clear();
    for (std::size_t i = 0; i < size; ++i)
    {
      std::pair<Key, T> pair;
      ar & pair;
      map.insert(std::move(pair));
    }
  }

  /**
   * @brief Serializes a flat map through split save/load functions.
   * @tparam Archive Archive type.
   * @tparam Key Map key type.
   * @tparam T Map mapped value type.
   * @tparam Compare Key comparison type.
   * @tparam Allocator Allocator type.
   * @param ar Archive used for serialization.
   * @param map Flat map to serialize.
   * @param version Serialization version.
   */
  template <class Archive, typename Key, typename T, typename Compare, typename Allocator>
  void serialize(
      Archive& ar, boost::container::flat_map<Key, T, Compare, Allocator>& map, const unsigned int version)
  {
    boost::serialization::split_free(ar, map, version);
  }
}

#endif // RODIN_SERIALIZATION_FLATMAP_H
