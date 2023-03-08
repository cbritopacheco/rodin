/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_DRAWABLES_LINE2D_H
#define RODIN_PLOT_BACKEND_DRAWABLES_LINE2D_H

#include <variant>
#include <algorithm>

#include <Eigen/Core>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/FlatGL.h>

#include "Rodin/Plot/Common.h"
#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Bases/BaseDrawable2D.h"


namespace Rodin::Plot::Backend::Renderer::Drawables
{
  class Line2D : public Bases::BaseDrawable2D
  {
   public:
   explicit Line2D(
      Object2D& object,
      DrawableGroup2D* group,
      const Eigen::ArrayX<float>& xData,
      const Eigen::ArrayX<float>& yData,
      const Magnum::Color3& color = Magnum::Color3{0, 0, 0},
      std::variant<LineStyle, DashTuple> lineStyle = Solid,
      float lineWidth = 5.0f,
      unsigned int lineSmoothness = 16
      );

   Line2D& setColor(const Magnum::Color4& color);
   Line2D& setLineStyle(const std::variant<LineStyle, DashTuple> lineStyle);
   Line2D& setLineWidth(float lineWidth);
   Line2D& setLineSmoothness(unsigned int lineSmoothness);

   private:
   void draw(
      const Magnum::Matrix3& transformationMatrix,
      Magnum::SceneGraph::Camera2D& camera) override;

   void computeMesh(Magnum::SceneGraph::Camera2D& camera);

   Eigen::ArrayX<float>                m_xData,
                                m_yData;
   Magnum::Color4                    m_color;
   std::variant<LineStyle,DashTuple>        m_lineStyle;
   float                          m_lineWidth;
   unsigned int                      m_lineSmoothness;

   // Renderer objects
   Magnum::GL::Mesh                   m_mesh;
   Magnum::SceneGraph::DrawableGroup2D      m_lines;
   Magnum::Shaders::FlatGL2D             m_flatShader;
  };
}

#endif
