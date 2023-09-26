/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_GUI_CURSOR_H
#define RODIN_PLOT_GUI_CURSOR_H

#include "Rodin/Plot/Common.h"

namespace Rodin::Plot::GUI
{
  class Cursor
  {
    public:
      enum SystemCursor
      {
        Arrow = SDL_SYSTEM_CURSOR_ARROW,
        IBeam = SDL_SYSTEM_CURSOR_IBEAM,
        Wait = SDL_SYSTEM_CURSOR_WAIT,
        Crosshair = SDL_SYSTEM_CURSOR_CROSSHAIR,
        WaitArrow = SDL_SYSTEM_CURSOR_WAITARROW,
        SizeNWSE = SDL_SYSTEM_CURSOR_SIZENWSE,
        SizeNESW = SDL_SYSTEM_CURSOR_SIZENESW,
        SizeWE = SDL_SYSTEM_CURSOR_SIZEWE,
        SizeAll = SDL_SYSTEM_CURSOR_SIZEALL,
        No = SDL_SYSTEM_CURSOR_NO,
        Hand = SDL_SYSTEM_CURSOR_HAND
      };

      Cursor() = default;

      Cursor(SystemCursor cursor)
      {
        m_cursor = SDL_CreateSystemCursor(static_cast<SDL_SystemCursor>(cursor));
      }

      Cursor(const Cursor&) = delete;
      void operator=(Cursor&) = delete;

      Cursor(Cursor&& other)
      {
        m_cursor = other.m_cursor;
        other.m_cursor = nullptr;
      }

      Cursor& operator=(Cursor&& other)
      {
        m_cursor = other.m_cursor;
        other.m_cursor = nullptr;
        return *this;
      }

      ~Cursor()
      {
        if (m_cursor)
          SDL_FreeCursor(m_cursor);
      }

      SDL_Cursor* getHandle()
      {
        return m_cursor;
      }

    private:
      SDL_Cursor* m_cursor;
  };
}
#endif
