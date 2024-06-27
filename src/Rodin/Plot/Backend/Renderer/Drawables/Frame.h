/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_DRAWABLES_FRAME_H
#define RODIN_PLOT_BACKEND_DRAWABLES_FRAME_H

#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Primitives/Square.h>
#include <Magnum/SceneGraph/Camera.h>

#include "Rodin/Plot/Backend/Bases/BaseDrawable2D.h"

namespace Rodin::Plot::Backend::Renderer::Drawables
{
  class Frame : public Bases::BaseDrawable2D
  {
    public:
      explicit
      Frame(
         Object2D& parent,
         DrawableGroup2D* group,
         const Magnum::Math::Vector<Scalar>2<float> bottomLeft,
         const Magnum::Math::Vector<Scalar>2<float> size,
         const Magnum::Color4& color = {0, 0, 0})
       : Bases::BaseDrawable2D(parent, group),
         m_color(color)
      {
        m_mesh = Magnum::MeshTools::compile(Magnum::Primitives::squareWireframe());
        getObject2D().scaleLocal({0.5, 0.5}).translate({0.5, 0.5});
        getObject2D().scale({size.x(), size.y()})
                     .translate({bottomLeft.x(), bottomLeft.y()});
      }

      Frame& setColor(const Magnum::Color4& color)
      {
        m_color = color;
        return *this;
      }

      Magnum::Color4 getColor() const
      {
        return m_color;
      }

    private:
      void draw(
        const Magnum::Matrix3& transformationMatrix,
        Magnum::SceneGraph::Camera2D& camera) override
      {
         Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::Multisampling);
         m_shader.setColor(m_color)
                 .setTransformationProjectionMatrix(
                    camera.projectionMatrix() * transformationMatrix)
                 .draw(m_mesh);
      }

      Magnum::Color4 m_color;
      Magnum::GL::Mesh m_mesh;
      Magnum::Shaders::FlatGL2D m_shader;
  };
}

#endif
