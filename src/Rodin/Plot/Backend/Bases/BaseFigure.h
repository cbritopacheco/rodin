/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef FEPLOT_BACKEND_BASES_BASEFIGURE_H
#define FEPLOT_BACKEND_BASES_BASEFIGURE_H

#include <string>
#include <memory>
#include <iostream>
#include <cassert>

#include <Magnum/Math/Vector2.h>

#include "Rodin/Core/Common.h"
#include "Rodin/Plot/Common.h"
#include "Rodin/Plot/Backend/Event.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseFigure
  {
   public:
    BaseFigure(const BaseFigure&) = delete;
    void operator=(const BaseFigure&) = delete;

    BaseFigure(const SDL_GLContext& glContext, const std::string& title,
       int width, int height);

    BaseFigure(BaseFigure&& other);

    virtual ~BaseFigure();

    void swapBuffers();

    void raise();

    FigureId getId() const;

    bool isVisible() const;

    std::string getTitle() const;

    BaseFigure& setVisible(bool isVisible);

    BaseFigure& setTitle(const std::string& title);

    bool makeCurrent();

    WindowHandle getWindowHandle();

    ConstWindowHandle getWindowHandle() const;

    Magnum::Math::Vector2<int> getWindowSize() const;

    Magnum::Math::Vector2<float> getDPIScaling() const;

    Magnum::Math::Vector2<int> getFrameBufferSize() const;

    virtual void drawContent() = 0;

    virtual void handle(const Backend::Event::MouseMotionEvent&) = 0;
    virtual void handle(const Backend::Event::MouseButtonEvent&) = 0;
    virtual void handle(const Backend::Event::MouseWheelEvent&) = 0;

   private:
    const SDL_GLContext&   m_glContext;
    WindowHandle        m_window;

    FigureId    m_id;
    std::string  m_title;
    int        m_width,
              m_height;
    bool       m_isVisible;

  };
}
#endif
