/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_BASES_BASECAMERA_H
#define RODIN_PLOT_BACKEND_BASES_BASECAMERA_H

#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Camera.h>

#include "Rodin/Plot/Backend/Renderer/Common.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseCamera3D
    : public Renderer::Object3D, public Magnum::SceneGraph::Camera3D
  {
    public:
      BaseCamera3D(Renderer::Object3D* parent)
        : Renderer::Object3D(parent),
          Magnum::SceneGraph::Camera3D(*this)
      {}
  };
}

#endif


