/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef PLOT_BACKEND_BASES_BASEEVENT_H
#define PLOT_BACKEND_BASES_BASEEVENT_H

#include <cstdint>

#include "Rodin/Plot/Common.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseEvent
  {
    public:
      BaseEvent(std::uint32_t timestamp) : m_timestamp(timestamp)
      {}

      std::uint32_t getTimestamp() const
      {
        return m_timestamp;
      }

      virtual ~BaseEvent()
      {}

    private:
      std::uint32_t m_timestamp;
  };
}

#endif
