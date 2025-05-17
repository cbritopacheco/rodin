#ifndef RODIN_SERIALIZATION_POLYTOPETYPE_H
#define RODIN_SERIALIZATION_POLYTOPETYPE_H

#include <type_traits>
#include <boost/serialization/array.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "Rodin/Array.h"

namespace boost::serialization
{
  template <class Archive, class ScalarType>
  void serialize(
      Archive & ar, const Rodin::Geometry::Polytope::Type& t, const unsigned int version)
  {
    ar & static_cast<std::underlying_type_t<Rodin::Geometry::Polytope::Type>>(t);
  }
}

#endif


