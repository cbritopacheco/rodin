#ifndef RODIN_SERIALIZATION_ARRAYSERIALIZATION_H
#define RODIN_SERIALIZATION_ARRAYSERIALIZATION_H

#include <boost/serialization/array.hpp>

#include "Rodin/Array.h"

namespace boost::serialization
{
  template <class Archive, class ScalarType>
  void serialize(
      Archive & ar, const Rodin::Array<ScalarType>& arr, const unsigned int version)
  {
    const size_t sz = arr.size();
    ar & sz;
    ar & boost::serialization::make_array(arr.data(), sz);
  }
}

#endif

