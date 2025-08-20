#ifndef RODIN_SERIALIZATION_FLATMAP_H
#define RODIN_SERIALIZATION_FLATMAP_H

#include <boost/container/flat_map.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/utility.hpp>

namespace boost::serialization
{
  template <class Archive, typename Key, typename T, typename Compare, typename Allocator>
  void save(
      Archive& ar, const boost::container::flat_map<Key, T, Compare, Allocator>& map, const unsigned int)
  {
    std::size_t size = map.size();
    ar & size;
    for (const auto& pair : map)
    {
      ar & pair;
    }
  }

  template <class Archive, typename Key, typename T, typename Compare, typename Allocator>
  void load(
      Archive& ar, boost::container::flat_map<Key, T, Compare, Allocator>& map, const unsigned int)
  {
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

  template <class Archive, typename Key, typename T, typename Compare, typename Allocator>
  void serialize(
      Archive& ar, boost::container::flat_map<Key, T, Compare, Allocator>& map, const unsigned int version)
  {
    boost::serialization::split_free(ar, map, version);
  }
}

#endif // RODIN_SERIALIZATION_FLATMAP_H
