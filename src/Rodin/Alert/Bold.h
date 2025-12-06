/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_COLOR_H
#define RODIN_ALERT_COLOR_H

namespace Rodin::Alert
{
  /**
   * @brief Empty tag type for bold text formatting.
   *
   * Tag type used to indicate bold text formatting in the Alert system.
   * This is typically used internally by the Text class attribute system.
   */
  struct BoldT {};

  /**
   * @brief Instance of BoldT tag type.
   *
   * Constant instance of the BoldT tag type for convenient usage.
   */
  static constexpr BoldT Bold;
}

#endif
