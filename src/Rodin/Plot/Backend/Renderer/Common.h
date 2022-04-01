/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_RENDERER_COMMON_H
#define RODIN_PLOT_BACKEND_RENDERER_COMMON_H

#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation2D.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>

namespace Rodin::Plot::Backend::Renderer
{
  using Object3D =
    Magnum::SceneGraph::Object<Magnum::SceneGraph::MatrixTransformation3D>;
  using Scene3D =
    Magnum::SceneGraph::Scene<Magnum::SceneGraph::MatrixTransformation3D>;
  using Object2D =
    Magnum::SceneGraph::Object<Magnum::SceneGraph::MatrixTransformation2D>;
  using Scene2D =
    Magnum::SceneGraph::Scene<Magnum::SceneGraph::MatrixTransformation2D>;
  using DrawableGroup2D = Magnum::SceneGraph::DrawableGroup2D;
  using DrawableGroup3D = Magnum::SceneGraph::DrawableGroup3D;
}

#endif
