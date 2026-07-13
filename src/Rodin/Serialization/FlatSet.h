/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FlatSet.h
 * @brief boost::serialization support for boost::container::flat_set.
 */
#ifndef RODIN_SERIALIZATION_FLATSET_H
#define RODIN_SERIALIZATION_FLATSET_H

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/nvp.hpp>

namespace boost::serialization
{
  /**
   * @brief Saves a flat set to an archive.
   * @tparam Archive Archive type.
   * @tparam T Set value type.
   * @tparam Compare Value comparison type.
   * @tparam Allocator Allocator type.
   * @param ar Archive used for serialization.
   * @param s Flat set to save.
   * @param version Serialization version.
   */
  template <class Archive, typename T, typename Compare, typename Allocator>
  void save(Archive& ar, const boost::container::flat_set<T, Compare, Allocator>& s,
    const unsigned int version)
  {
    (void)version;
    std::size_t sz = s.size();
    ar& BOOST_SERIALIZATION_NVP(sz);
#if BOOST_VERSION >= 108100
     if(sz > 0)
     {
        const auto & seq = s.get_sequence_cref();
        ar.save_binary(reinterpret_cast<const char*>(seq.data()), sz * sizeof(T));
     }
#else
     for(const auto& elem : s)
        ar & BOOST_SERIALIZATION_NVP(elem);
#endif
  }

  /**
   * @brief Loads a flat set from an archive.
   * @tparam Archive Archive type.
   * @tparam T Set value type.
   * @tparam Compare Value comparison type.
   * @tparam Allocator Allocator type.
   * @param ar Archive used for serialization.
   * @param s Flat set to load.
   * @param version Serialization version.
   */
  template <class Archive, typename T, typename Compare, typename Allocator>
  void load(Archive& ar, boost::container::flat_set<T, Compare, Allocator>& s,
    const unsigned int version)
  {
    (void)version;
    std::size_t sz;
    ar& BOOST_SERIALIZATION_NVP(sz);
    s.clear();
#if BOOST_VERSION >= 108100
     if(sz > 0)
     {
        auto & seq = s.get_sequence_ref();
        seq.resize(sz);
        ar.load_binary(reinterpret_cast<char*>(seq.data()),
                       sz * sizeof(T));
     }
#else
     for(std::size_t i = 0; i < sz; ++i)
     {
        T elem;
        ar & BOOST_SERIALIZATION_NVP(elem);
        s.insert(elem);
     }
#endif
  }

  /**
   * @brief Serializes a flat set through split save/load functions.
   * @tparam Archive Archive type.
   * @tparam T Set value type.
   * @tparam Compare Value comparison type.
   * @tparam Allocator Allocator type.
   * @param ar Archive used for serialization.
   * @param s Flat set to serialize.
   * @param version Serialization version.
   */
  template <class Archive, typename T, typename Compare, typename Allocator>
  void serialize(Archive & ar,
                 boost::container::flat_set<T, Compare, Allocator>& s,
                 const unsigned int version)
  {
    boost::serialization::split_free(ar, s, version);
  }
}

#endif
