/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEEVENT_H
#define PLOT_BACKEND_EVENT_MOUSEEVENT_MOUSEEVENT_H

#include "Rodin/Plot/Backend/Bases/BaseEvent.h"

namespace Rodin::Plot::Backend::Event
{
  class MouseEvent : public Backend::Bases::BaseEvent
  {
    public:
    enum Button
    {
      LEFT    = SDL_BUTTON_LEFT,
      MIDDLE  = SDL_BUTTON_MIDDLE,
      RIGHT   = SDL_BUTTON_RIGHT,
      X1      = SDL_BUTTON_X1,
      X2      = SDL_BUTTON_X2
    };

    MouseEvent(std::uint32_t timestamp)
      : BaseEvent(timestamp)
    {}
  };
}

#endif
