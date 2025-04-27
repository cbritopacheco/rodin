/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <string>

#include <Magnum/GL/DefaultFramebuffer.h>

#include "Rodin/Plot/Plot.h"

#include "Figure.h"

namespace Rodin::Plot::Artist
{
  unsigned int Figure::s_figureCount = 1;

  Figure::Figure(Plot& pltRef)
    : Figure(pltRef, 640, 480)
  {}

  Figure::Figure(Plot& pltRef, int width, int height)
    : Figure(pltRef, width, height, std::string("Figure ") +
        std::to_string(s_figureCount++))
  {}

  Figure::Figure(Plot& pltRef,
      int width, int height, const std::string& title)
    : BaseFigure(pltRef.getGLContext(), title, width, height),
    m_pltRef(pltRef),
    m_camera(m_scene.addChild<Backend::Renderer::Object2D>()),
    m_redraw(true)
  {
    // Set size of view space equal to the figure space
    m_camera.setProjectionMatrix(
        Magnum::Math::Matrix3<float>::projection({
          static_cast<float>(width),
          static_cast<float>(height)
          })
        );

    // Origin is the lower left corner
    m_camera.getObject2D().translate({width / 2.0f, height / 2.0f});
  }

  void Figure::drawContent()
  {
    for (auto& ax : m_axes)
    {
      auto windowSize = getWindowSize();
      auto frameBufferSize = getFrameBufferSize();
      auto xRatio = frameBufferSize.x() / static_cast<float>(windowSize.x()),
           yRatio = frameBufferSize.y() / static_cast<float>(windowSize.y());

      Magnum::GL::defaultFramebuffer.setViewport({
          {
          static_cast<int>(xRatio * ax.axes->getBottomLeft().x()),
          static_cast<int>(yRatio * ax.axes->getBottomLeft().y())
          },
          {
          static_cast<int>(xRatio * (
                ax.axes->getBottomLeft().x() + ax.axes->getSize().x())),
          static_cast<int>(yRatio * (
                ax.axes->getBottomLeft().y() + ax.axes->getSize().y()))
          }
          });

      ax.axes->drawContent();

      m_redraw = false;
    }

    // Reset the viewport to draw the Drawables which make up the graphical
    // interface of the window
    auto size = getFrameBufferSize();
    Magnum::GL::defaultFramebuffer.setViewport({
        {0, 0},
        {size.x(), size.y()}
        });
    m_camera.draw(getDrawableGroup());
  }

  void Figure::handle(const Backend::Event::MouseMotionEvent& e)
  {
    m_redraw = true;
    for (auto it = m_axes.rbegin(); it != m_axes.rend(); it++)
    {
      auto& ax = *it;
      if (ax.axes->getBoundingBox().contains({e.getX(), e.getY()}))
      {
        getPlot().setCursor(GUI::Cursor::Crosshair);
        ax.axes->handle(e);
        return;
      }
    }
    getPlot().setCursor(GUI::Cursor::Arrow);
  }

  void Figure::handle(const Backend::Event::MouseButtonEvent& e)
  {
    m_redraw = true;
    for (auto it = m_axes.rbegin(); it != m_axes.rend(); it++)
    {
      auto& ax = *it;
      if (ax.axes->getBoundingBox().contains({e.getX(), e.getY()}))
      {
        switch (e.getButtonState())
        {
          case Backend::Event::MouseButtonEvent::PRESSED:
          getPlot().setCursor(GUI::Cursor::SizeAll);
          break;
          case Backend::Event::MouseButtonEvent::RELEASED:
          getPlot().setCursor(GUI::Cursor::Crosshair);
          break;
        }
        ax.axes->handle(e);
        return;
      }
    }
    getPlot().setCursor(GUI::Cursor::Arrow);
  }

  void Figure::handle(const Backend::Event::MouseWheelEvent& e)
  {
    m_redraw = true;
    auto mousePos = getPlot().getMousePosition();
    for (auto it = m_axes.rbegin(); it != m_axes.rend(); it++)
    {
      auto& ax = *it;
      if (ax.axes->getBoundingBox().contains(mousePos.x(), mousePos.y()))
      {
        ax.axes->handle(e);
        return;
      }
    }
  }

  void Figure::clear()
  {
    Magnum::GL::Renderer::setClearColor(Magnum::Color3{1.0f, 1.0f, 1.0f});
    Magnum::GL::defaultFramebuffer.clear(Magnum::GL::FramebufferClear::Color |
        Magnum::GL::FramebufferClear::Depth);
  }


  Plot& Figure::getPlot()
  {
    return m_pltRef;
  }

  const Plot& Figure::getPlot() const
  {
    return m_pltRef;
  }

  Backend::Renderer::Object2D& Figure::getObject2D()
  {
    return m_scene;
  }

  Backend::Renderer::DrawableGroup2D& Figure::getDrawableGroup()
  {
    return m_drawables;
  }

  const Backend::Renderer::DrawableGroup2D& Figure::getDrawableGroup() const
  {
    return m_drawables;
  }

  bool Figure::redraw() const
  {
    return m_redraw;
  }
}
