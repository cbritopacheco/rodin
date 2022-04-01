/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "BaseFigure.h"

namespace Rodin::Plot::Backend::Bases
{
  BaseFigure::BaseFigure(
      const SDL_GLContext& glContext,
      const std::string& title,
      int width,
      int height)
    : m_glContext(glContext),
      m_title(title), m_width(width),
      m_height(height), m_isVisible(false)
  {
    m_window = SDL_CreateWindow(
        title.c_str(),
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        width,
        height,
        SDL_WINDOW_OPENGL | SDL_WINDOW_HIDDEN | SDL_WINDOW_ALLOW_HIGHDPI
        );
      if (m_window == NULL)
      {
        printf("Could not create window: %s\n", SDL_GetError());
        exit(EXIT_FAILURE);
      } else {
        m_id = SDL_GetWindowID(m_window);
      }
  }

  BaseFigure::BaseFigure(BaseFigure&& other)
    : m_glContext(other.m_glContext),
      m_window(other.m_window),
      m_id(other.m_id),
      m_title(std::move(other.m_title)),
      m_width(other.m_width),
      m_height(other.m_height),
      m_isVisible(other.m_isVisible)
  {
    other.m_window = nullptr;
    other.m_id = 0;
    other.m_title = "";
    other.m_width = 0;
    other.m_height = 0;
    other.m_isVisible = false;
  }

  BaseFigure::~BaseFigure()
  {
    if (m_window)
      SDL_DestroyWindow(m_window);
  }

  void BaseFigure::swapBuffers()
  {
    SDL_GL_SwapWindow(m_window);
  }

  void BaseFigure::raise()
  {
    SDL_RaiseWindow(m_window);
  }

  FigureId BaseFigure::getId() const
  {
    return m_id;
  }

  bool BaseFigure::isVisible() const
  {
    return m_isVisible;
  }

  BaseFigure& BaseFigure::setVisible(bool isVisible)
  {
    m_isVisible = isVisible;
    if (isVisible)
      SDL_ShowWindow(m_window);
    else
      SDL_HideWindow(m_window);
    return *this;
  }

  std::string BaseFigure::getTitle() const
  {
    return m_title;
  }

  BaseFigure& BaseFigure::setTitle(const std::string& title)
  {
    SDL_SetWindowTitle(m_window, title.c_str());
    m_title = title;
    return *this;
  }

  WindowHandle BaseFigure::getWindowHandle()
  {
    return m_window;
  }

  ConstWindowHandle BaseFigure::getWindowHandle() const
  {
    return m_window;
  }

  bool BaseFigure::makeCurrent()
  {
    if (SDL_GL_MakeCurrent(m_window, m_glContext) < 0)
    {
      std::cerr << "Failed to make the context current" << std::endl;
      return false;
    }
    else
      return true;
  }

  Eigen::Array2<int> BaseFigure::getFrameBufferSize() const
  {
    Eigen::Array2<int> size;
    SDL_GL_GetDrawableSize(m_window, &size.x(), &size.y());
    return size;
  }

  Eigen::Array2<int> BaseFigure::getWindowSize() const
  {
    Eigen::Array2<int> size;
    SDL_GetWindowSize(m_window, &size.x(), &size.y());
    return size;
  }

  Eigen::Array2<float> BaseFigure::getDPIScaling() const
  {
    Eigen::Array2<float> dpi;
    if(SDL_GetDisplayDPI(0, nullptr, &dpi.x(), &dpi.y()) == 0)
    {
      return dpi;
    }
    else
    {
      std::cerr << "Failed to make the context current" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}
