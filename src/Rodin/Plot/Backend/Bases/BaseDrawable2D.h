/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_BASES_BASEDRAWABLE2D_H
#define RODIN_PLOT_BACKEND_BASES_BASEDRAWABLE2D_H

#include <Magnum/SceneGraph/Drawable.h>

#include "Rodin/Plot/Backend/Renderer/Common.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseDrawable2D : public Magnum::SceneGraph::Drawable2D
  {
   public:
    BaseDrawable2D(
       Renderer::Object2D& object, Renderer::DrawableGroup2D* group)
    : Magnum::SceneGraph::Drawable2D{object, group},
      m_object(object)
   {}

   Renderer::Object2D& getObject2D()
   {
    return m_object;
   }

   const Renderer::Object2D& getObject2D() const
   {
    return m_object;
   }

   private:
    Renderer::Object2D& m_object;
  };
}

#endif
