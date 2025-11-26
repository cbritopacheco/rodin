/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_CAST_H
#define RODIN_UTILITY_CAST_H

/**
 * @file
 * @brief Defines the Cast template for type conversions.
 */

namespace Rodin::Utility
{
  /**
   * @brief Template for defining type conversion operations.
   * @ingroup UtilityModule
   * @tparam From The source type to convert from.
   * @tparam To The destination type to convert to.
   *
   * The Cast template provides a customization point for defining
   * explicit type conversions between arbitrary types. Specializations
   * of this template can implement custom conversion logic for specific
   * type pairs.
   *
   * This is useful when standard conversions are not available or when
   * custom conversion semantics are required.
   */
  template <class From, class To>
  class Cast;
}

#endif
