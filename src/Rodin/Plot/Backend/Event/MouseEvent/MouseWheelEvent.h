/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEWHEELEVENT_H
#define PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEWHEELEVENT_H

#include "MouseEvent.h"

namespace Rodin::Plot::Backend::Event
{
  class MouseWheelEvent : public MouseEvent
  {
   public:
    enum WheelDirection
    {
      NORMAL = SDL_MOUSEWHEEL_NORMAL,
      FLIPPED = SDL_MOUSEWHEEL_FLIPPED
    };

   MouseWheelEvent(std::uint32_t timestamp, std::int32_t x, std::int32_t y,
      WheelDirection wheelDirection)
    : MouseEvent(timestamp), m_x(x), m_y(y), m_wheelDirection(wheelDirection)
   {}

   std::int32_t getX() const
   {
    return m_x;
   }

   std::int32_t getY() const
   {
    return m_y;
   }

   WheelDirection getWheelDirection() const
   {
    return m_wheelDirection;
   }

   private:
   std::int32_t    m_x,
               m_y;
   WheelDirection   m_wheelDirection;
  };
}

#endif
