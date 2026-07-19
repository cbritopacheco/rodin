/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file UnorderedMap.h
 * @brief boost::serialization support for boost::unordered_map.
 */
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
