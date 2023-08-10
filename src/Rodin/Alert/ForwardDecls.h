/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_FORWARDDECLS_H
#define RODIN_ALERT_FORWARDDECLS_H

#include <cstdlib>

namespace Rodin::Alert
{
  template <class Prefix>
  class Message;

  class ResetT;

  class StylizeT;

  template <size_t RED, size_t GREEN, size_t BLUE>
  struct RGB;

  template <class Code>
  class Color;
}

#endif

