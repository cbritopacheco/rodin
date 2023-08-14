/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_BASES_BASECAMERA2D_H
#define RODIN_PLOT_BACKEND_BASES_BASECAMERA2D_H

#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Camera.h>

#include "Rodin/Plot/Backend/Renderer/Common.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseCamera2D : public Magnum::SceneGraph::Camera2D
  {
   public:
    BaseCamera2D(Renderer::Object2D& object)
      : Magnum::SceneGraph::Camera2D(object),
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


