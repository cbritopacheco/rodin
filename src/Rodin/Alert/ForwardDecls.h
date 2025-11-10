/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_FORWARDDECLS_H
#define RODIN_ALERT_FORWARDDECLS_H

#include <cstdlib>

/**
 * @file ForwardDecls.h
 * @brief Forward declarations for Alert module types.
 * @ingroup AlertModule
 *
 * This header contains forward declarations of key types used throughout
 * the Alert module to minimize header dependencies and improve compilation
 * times.
 */

namespace Rodin::Alert
{
  /**
   * @brief Forward declaration of the Message template class.
   * @tparam Prefix The prefix type for the message.
   *
   * Base template class for formatted message output with customizable prefixes.
   */
  template <class Prefix>
  class Message;

  /**
   * @brief Forward declaration of ResetT tag type.
   *
   * Tag type for resetting terminal text formatting and colors.
   */
  struct ResetT;

  /**
   * @brief Forward declaration of StylizeT tag type.
   *
   * Tag type for enabling terminal color output.
   */
  struct StylizeT;

  /**
   * @brief Forward declaration of RGB color template.
   * @tparam RED Red color component (0-255).
   * @tparam GREEN Green color component (0-255).
   * @tparam BLUE Blue color component (0-255).
   *
   * Template struct representing RGB color values for custom terminal colors.
   */
  template <size_t RED, size_t GREEN, size_t BLUE>
  struct RGB;

  /**
   * @brief Forward declaration of Color template class.
   * @tparam Code The color code type.
   *
   * Template class for custom terminal colors with type-safe color codes.
   */
  template <class Code>
  class Color;
}

#endif

