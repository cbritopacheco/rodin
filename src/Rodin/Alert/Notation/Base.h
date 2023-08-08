/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_NOTATION_BASE_H
#define RODIN_ALERT_NOTATION_BASE_H

#include "Rodin/Alert/ForwardDecls.h"

namespace Rodin::Alert::Notation
{
  class Base
  {
    public:
      virtual const char* toCharString() const = 0;
  };

  class ArrowT final : public Base
  {
    public:
      const char* toCharString() const override
      {
        return s_charString;
      }

    private:
      static const char* s_charString;
  };

  static constexpr ArrowT Arrow;
}

#endif
