#ifndef RODIN_SERIALIZATION_ARRAYSERIALIZATION_H
#define RODIN_SERIALIZATION_ARRAYSERIALIZATION_H

#include <boost/serialization/array.hpp>
#include <boost/serialization/boost_array.hpp>

#include "Rodin/Array.h"

namespace boost::serialization
{
  template <class Archive, typename ScalarType>
  void save(Archive & ar,
            const Rodin::Array<ScalarType>& arr,
            const unsigned int)
  {
    const size_t sz = arr.size();
    ar & sz;
    ar & boost::serialization::make_array(arr.data(), sz);
  }

  template <class Archive, typename ScalarType>
  void load(Archive & ar,
            Rodin::Array<ScalarType>& arr,
            const unsigned int)
  {
    size_t sz;
    ar & sz;
    arr.resize(sz);
    ar & boost::serialization::make_array(arr.data(), sz);
  }

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

