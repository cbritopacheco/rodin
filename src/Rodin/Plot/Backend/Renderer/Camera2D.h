/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_RENDERER_CAMERA2D_H
#define RODIN_PLOT_BACKEND_RENDERER_CAMERA2D_H

#include "Rodin/Plot/Backend/Renderer/Common.h"
#include "Rodin/Plot/Backend/Bases/BaseCamera2D.h"

namespace Rodin::Plot::Backend::Renderer
{
  class Camera2D : public Bases::BaseCamera2D
  {
    public:
      Camera2D(Object2D& scene)
        : Bases::BaseCamera2D(scene)
      {}
  };
}

#endif
