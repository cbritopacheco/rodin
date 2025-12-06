/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CONTEXT_BASE_H
#define RODIN_CONTEXT_BASE_H

#include <boost/serialization/access.hpp>

#include "Rodin/Configure.h"

namespace Rodin::Context
{
  /**
   * @brief Abstract base class for all execution contexts.
   *
   * Base provides the fundamental interface that all execution contexts
   * must implement. It defines the basic serialization interface required
   * for context management and provides a foundation for context-specific
   * optimizations and resource management.
   */
  class Base
  {
    friend class boost::serialization::access;

    public:
      /**
       * @brief Serialization method for Boost.Serialization support.
       *
       * @tparam Archive Archive type for serialization
       * @param ar Archive instance
       * @param version Serialization version number
       *
       * Provides serialization support for context objects, enabling
       * context state to be saved and restored across different execution
       * sessions or distributed across computing nodes.
       */
      template <class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {}
  };
}

#endif


