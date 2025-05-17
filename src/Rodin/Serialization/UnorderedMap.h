#ifndef RODIN_SERIALIZATION_UNORDEREDMAP_H
#define RODIN_SERIALIZATION_UNORDEREDMAP_H

#include <boost/unordered_map.hpp>
#include <boost/serialization/unordered_map.hpp>

#if BOOST_VERSION < 108400
#include <boost/serialization/boost_unordered_map.hpp>
#endif

#endif

