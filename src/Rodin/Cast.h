/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CAST_H
#define RODIN_CAST_H

/**
 * @file
 * @brief Defines the Cast template for type conversion operations.
 */

namespace Rodin
{
  /**
   * @defgroup RodinCasting Type Casting Utilities
   * @brief Safe and efficient type casting operations.
   * 
   * This module provides utilities for safe type conversions between
   * related types in the Rodin framework. The casting operations are
   * designed to maintain type safety while providing efficient conversions
   * for polymorphic objects and template instantiations.
   */

  /**
   * @ingroup RodinCasting
   * @brief Template class for type casting operations.
   *
   * Cast provides a generic interface for converting between types,
   * with specializations for specific type pairs. This enables safe
   * and efficient type conversions while maintaining compile-time
   * type checking.
   *
   * @tparam From Source type to cast from
   * @tparam To Target type to cast to
   *
   * ## Usage
   * ```cpp
   * // Example casting operation
   * Cast<SourceType, TargetType> converter;
   * TargetType result = converter(sourceObject);
   * ```
   */
  template <class From, class To>
  class Cast;
}

#endif

