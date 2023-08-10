/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_RAISE_H
#define RODIN_ALERT_RAISE_H

namespace Rodin::Alert
{
  /**
   * @brief Empty class tag type to raise an Alert.
   *
   * RaiseT is an empty class tag type to specify that a derived object of
   * Alert should be raised.
   */
  struct RaiseT
  {
    explicit constexpr RaiseT() = default;
  };

  /**
   * @brief Instance of RaiseT
   *
   * Instance of the empty struct tag type RaiseT.
   */
  static constexpr RaiseT Raise;
}

#endif

