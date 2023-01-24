/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_FORWARDDECLS_H
#define RODIN_PLOT_FORWARDDECLS_H

namespace Rodin::Plot
{
  namespace Backend
  {
   namespace Bases
   {
    class BaseArtist2D;
    class BaseTopLevelArtist2D;
    class BaseAxes;
    class BaseFigure;
    class BaseEvent;
    class BaseDrawable2D;
   }

   namespace Renderer
   {
    namespace Drawables
    {
      class Line2D;
    }
   }
  }

  class Plot;

  namespace Artist
  {
   class Figure;

   namespace Axes
   {
    class Axes2D;
   }

   namespace Lines
   {
    class Line2D;
   }
  }
}

#endif
