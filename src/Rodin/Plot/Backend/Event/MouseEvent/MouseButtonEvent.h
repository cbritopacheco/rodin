/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEBUTTONEVENT_H
#define PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEBUTTONEVENT_H

#include "MouseEvent.h"

namespace Rodin::Plot::Backend::Event
{
  class MouseButtonEvent : public MouseEvent
  {
   public:
    enum ButtonState
    {
      PRESSED = SDL_PRESSED,
      RELEASED = SDL_RELEASED
    };

    enum Clicks
    {
      SINGLE_CLICK = 1,
      DOUBLE_CLICK = 2
    };

    MouseButtonEvent(
       std::uint32_t timestamp,
       Button button,
       ButtonState buttonState,
       Clicks clicks,
       std::int32_t x,
       std::int32_t y
       )
      :  MouseEvent(timestamp),
        m_button(button),
        m_buttonState(buttonState),
        m_clicks(clicks),
        m_x(x),
        m_y(y)
    {}

    Button getButton() const
    {
      return m_button;
    }

    ButtonState getButtonState() const
    {
      return m_buttonState;
    }

    Clicks getClicks() const
    {
      return m_clicks;
    }

    std::int32_t getX() const
    {
      return m_x;
    }

    std::int32_t getY() const
    {
      return m_y;
    }

   private:
    Button      m_button;
    ButtonState  m_buttonState;
    Clicks      m_clicks;
    std::int32_t  m_x,
              m_y;
  };
}
#endif
