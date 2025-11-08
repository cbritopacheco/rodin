#ifndef RODIN_SERIALIZATION_UNORDEREDMAP_H
#define RODIN_SERIALIZATION_UNORDEREDMAP_H

// Must come first so library_version_type is visible to Boost helpers
#include <boost/serialization/library_version_type.hpp>
#include <boost/serialization/version.hpp>

// Legacy containers: boost::unordered_map (namespace boost)
#include <boost/unordered_map.hpp>
#include <boost/serialization/unordered_map.hpp>

// New containers: boost::unordered::unordered_map (namespace boost::unordered)
#include <boost/unordered/unordered_map.hpp>

#if BOOST_VERSION < 108400
  #include <boost/serialization/boost_unordered_map.hpp>
#endif

#endif // RODIN_SERIALIZATION_UNORDEREDMAP_H
