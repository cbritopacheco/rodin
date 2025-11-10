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
   * @brief Empty tag type for raising an Alert message.
   * @ingroup AlertModule
   *
   * Tag type used to trigger the raising (output) of an Alert message.
   * When streamed into a Message object using operator<<, it causes
   * the message to be displayed and potentially triggers other actions
   * (such as program termination for Exception types).
   *
   * Example usage:
   * @code
   * Exception() << "An error occurred" << Raise;
   * @endcode
   */
  struct RaiseT
  {
    explicit constexpr RaiseT() = default;
  };

  /**
   * @brief Instance of RaiseT tag type.
   * @ingroup AlertModule
   *
   * Constant instance of the RaiseT tag type for convenient usage.
   * Stream this into a Message to trigger its raise() method.
   */
  static constexpr RaiseT Raise;
}

#endif

