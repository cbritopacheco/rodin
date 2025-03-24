#ifndef RODIN_SERIALIZATION_FLATSET_H
#define RODIN_SERIALIZATION_FLATSET_H

#include <boost/serialization/split_free.hpp>

#include "Rodin/FlatSet.h"

namespace boost::serialization
{
  template <class Archive, typename T, typename Compare, typename Allocator>
  void save(Archive & ar,
            const boost::container::flat_set<T, Compare, Allocator>& s,
            const unsigned int)
  {
    size_t sz = s.size();
    ar & sz;
    if(sz > 0)
    {
      const auto& sequence = s.get_sequence_cref();
      ar.save_binary(
          reinterpret_cast<const char*>(sequence.data()), sz * sizeof(T));
    }
  }

  template <class Archive, typename T, typename Compare, typename Allocator>
  void load(Archive & ar,
            boost::container::flat_set<T, Compare, Allocator>& s,
            const unsigned int)
  {
    size_t sz = s.size();
    ar & sz;
    s.clear();
    if(sz > 0)
    {
      auto& sequence = s.get_sequence_ref();
      sequence.resize(sz);
      ar.load_binary(
          reinterpret_cast<char*>(sequence.data()), sz * sizeof(T));
    }
  }

  template <class Archive, typename T, typename Compare, typename Allocator>
  void serialize(
      Archive & ar, boost::container::flat_set<T, Compare, Allocator>& s, const unsigned int version)
  {
    boost::serialization::split_free(ar, s, version);
  }
}

#endif

