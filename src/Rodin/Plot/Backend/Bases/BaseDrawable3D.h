/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_BASES_BASEDRAWABLE3D_H
#define RODIN_PLOT_BACKEND_BASES_BASEDRAWABLE3D_H

#include <Magnum/SceneGraph/Drawable.h>

#include "Rodin/Plot/Backend/Renderer/Common.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseDrawable3D
    : public Renderer::Object3D, public Magnum::SceneGraph::Drawable3D
  {
    public:
      BaseDrawable3D(
          Renderer::Object3D* parent,
          Magnum::SceneGraph::DrawableGroup3D* group
          )
      : Renderer::Object3D{parent},
        Magnum::SceneGraph::Drawable3D{*this, group}
    {}
  };
}

#endif
