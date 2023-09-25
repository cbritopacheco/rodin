/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEMOTIONEVENT_H
#define PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEMOTIONEVENT_H

#include "MouseEvent.h"

namespace Rodin::Plot::Backend::Event
{
  class MouseMotionEvent : public MouseEvent
  {
    public:
      MouseMotionEvent(std::uint32_t timestamp, std::uint32_t buttonState,
          std::int32_t x, std::int32_t y, std::int32_t motionX,
          std::int32_t motionY)
        : MouseEvent(timestamp),
        m_buttonState(buttonState), m_x(x), m_y(y),
        m_motionX(motionX), m_motionY(motionY)
      {}

      bool isButtonPressed(Button button) const
      {
        std::uint32_t state;
        switch(button)
        {
          case LEFT:
          state = SDL_BUTTON_LMASK;
          break;
          case MIDDLE:
          state = SDL_BUTTON_MMASK;
          break;
          case RIGHT:
          state = SDL_BUTTON_LMASK;
          break;
          case X1:
          state = SDL_BUTTON_X1MASK;
          break;
          case X2:
          state = SDL_BUTTON_X2MASK;
          break;
        }
        return m_buttonState & state;
      }

      std::int32_t getX() const
      {
        return m_x;
      }

      std::int32_t getY() const
      {
        return m_y;
      }

      std::int32_t getMotionX() const
      {
        return m_motionX;
      }

      std::int32_t getMotionY() const
      {
        return m_motionY;
      }

    private:
      std::uint32_t m_buttonState;
      std::int32_t  m_x,
                    m_y,
                    m_motionX,
                    m_motionY;
  };
}

#endif
