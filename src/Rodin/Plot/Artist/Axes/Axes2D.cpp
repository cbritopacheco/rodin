/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/Math/Matrix3.h>

#include "Rodin/Plot/Plot.h"
#include "Rodin/Plot/Backend/Renderer/Drawables/Line2D.h"

#include "Axes2D.h"

namespace Rodin::Plot::Artist::Axes
{
  Axes2D::Axes2D(
      Figure& fig,
      Magnum::Math::Vector2<int> bottomLeft,
      Magnum::Math::Vector2<int> size,
      bool frameEnabled)
    : BaseAxes(fig, bottomLeft, size, frameEnabled),
      m_camera(m_scene.addChild<Backend::Renderer::Object2D>())
  {
    setXLimits({-1, 1});
    setYLimits({-1, 1});
  }

  void Axes2D::drawContent()
  {
    m_camera.setProjectionMatrix(
        Magnum::Math::Matrix3<float>::projection({
            getXLimits().right - getXLimits().left,
            getYLimits().top - getYLimits().bottom
          })
        );
    m_camera.getObject2D().setTransformation(
        Magnum::Math::Matrix3<float>::translation(
        {
          (getXLimits().right + getXLimits().left) / 2.0f,
          (getYLimits().top + getYLimits().bottom) / 2.0f
        })
      );
    m_camera.draw(getDrawableGroup());
  }

  void Axes2D::handle(const Backend::Event::MouseMotionEvent& e)
  {
    if (e.isButtonPressed(Backend::Event::MouseEvent::LEFT))
    {
      auto xlim = getXLimits();
      auto ylim = getYLimits();
      auto axesSize = getSize();

      float dx = e.getMotionX() * (xlim.right - xlim.left) / axesSize.x(),
            dy = e.getMotionY() * (ylim.top - ylim.bottom) / axesSize.y();

      xlim.left -= dx;
      xlim.right -= dx;
      ylim.bottom -= dy;
      ylim.top -= dy;

      setXLimits(xlim);
      setYLimits(ylim);
      getFigure().getPlot().setCursor(GUI::Cursor::SizeAll);
    }
  }

  void Axes2D::handle(const Backend::Event::MouseButtonEvent& e)
  {
  }

  void Axes2D::handle(const Backend::Event::MouseWheelEvent& e)
  {
    auto xlim = getXLimits();
    auto ylim = getYLimits();
    auto axesSize = getSize();

    float dx = e.getY() * (xlim.right - xlim.left) / axesSize.x(),
          dy = e.getY() * (ylim.top - ylim.bottom) / axesSize.y();

    xlim.left   -= dx;
    xlim.right  += dx;
    ylim.bottom -= dy;
    ylim.top    += dy;

    setXLimits(xlim);
    setYLimits(ylim);
  }

  Backend::Renderer::Object2D& Axes2D::getObject2D()
  {
    return m_scene;
  }

  Backend::Renderer::DrawableGroup2D& Axes2D::getDrawableGroup()
  {
    return m_drawables;
  }

  const Backend::Renderer::DrawableGroup2D& Axes2D::getDrawableGroup() const
  {
    return m_drawables;
  }
}
