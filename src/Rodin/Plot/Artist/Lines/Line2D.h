/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_ARTIST_LINES_LINE2D_H
#define RODIN_PLOT_ARTIST_LINES_LINE2D_H

#include "Rodin/Plot/Common.h"
#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Bases/BaseArtist2D.h"
#include "Rodin/Plot/Backend/Renderer/Drawables/Line2D.h"

namespace Rodin::Plot::Artist::Lines
{
  class Line2D : public Backend::Bases::BaseArtist2D
  {
   public:
   template <class ... Args>
   Line2D(Backend::Bases::BaseArtist2D& parent, Args&&... args)
    : Backend::Bases::BaseArtist2D(parent),
      m_object(parent.getObject2D().addChild<Backend::Renderer::Object2D>())
   {
    draw<Backend::Renderer::Drawables::Line2D>(std::forward<Args>(args)...);
   }

   virtual Backend::Renderer::Object2D& getObject2D() override
   {
    return m_object;
   }

   private:
   Backend::Renderer::Object2D&        m_object;
  };
}

#endif
