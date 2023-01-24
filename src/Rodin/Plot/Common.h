/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_COMMON_H
#define RODIN_PLOT_COMMON_H

#include <vector>

#include <SDL.h>

namespace Rodin::Plot
{
  /**
  * Represents the unique identifier of an Artist::Figure instance.
  */
  using FigureId       = decltype(SDL_GetWindowID(std::declval<SDL_Window*>()));

  using WindowHandle    = SDL_Window*;
  using ConstWindowHandle = const SDL_Window*;

  enum LineStyle
  {
   Solid,
   Dashed,
   DashDotted,
   Dotted
  };

  class DashTuple
  {
   public:
   DashTuple(int offset, std::vector<int> sequence)
    : m_offset(offset), m_sequence(sequence)
   {}

   const std::vector<int>& getSequence() const
   {
    return m_sequence;
   }

   int getOffset()
   {
    return m_offset;
   }

   private:
    int m_offset;
    std::vector<int> m_sequence;
  };
}


#endif
