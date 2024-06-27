/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Magnum/GL/Renderer.h>

#include "Rodin/Plot/Backend/Event.h"

#include "Plot.h"

namespace Rodin::Plot
{
  Plot::Plot()
    : Plot(Configuration{})
    {}

  Plot::Plot(const Configuration& configuration)
    : m_minimalLoopPeriod(0),
      m_isVsyncEnabled(true)
  {
    //Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
      SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
      std::exit(EXIT_FAILURE);
    }
    assert(SDL_WasInit(SDL_INIT_VIDEO) == SDL_INIT_VIDEO);

    // TODO: OpenGL versioning
    SDL_GL_SetAttribute( SDL_GL_CONTEXT_PROFILE_MASK,
        SDL_GL_CONTEXT_PROFILE_CORE );

    // Enable double buffering
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    // Enable anti-aliasing
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS,
        configuration.getSampleCount() > 1 ? 1 : 0);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, configuration.getSampleCount());

    // Now we will create the underlying SDL OpenGL context and the Magnum
    // context. All figures/windows will share the same context.
    m_context.emplace(Corrade::NoCreate, 0, nullptr);

    // We need a dummy window to create the OpenGL context. This should not
    // pose a problem since the OpenGL context should be independent of the
    // drawable.
    m_initWindow = SDL_CreateWindow(
        nullptr, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        1, 1, SDL_WINDOW_OPENGL | SDL_WINDOW_HIDDEN);

    // Create underlying context
    m_glContext = SDL_GL_CreateContext(m_initWindow);
    if (m_glContext == NULL)
    {
      printf("Context creation failed: %s", SDL_GetError());
      std::exit(EXIT_FAILURE);
    }
    else
    {
      // Make the context current for the initial window and create the Magnum
      // OpenGL context
      SDL_GL_MakeCurrent(m_initWindow, m_glContext);
      Magnum::GL::Context::makeCurrent(nullptr);
      if(!m_context->tryCreate())
      {
        printf("Magnum context creation failed!\n");
        std::exit(EXIT_FAILURE);
      }
      Magnum::GL::Context::makeCurrent(&*m_context);
    }

    m_arrowCursor = GUI::Cursor(GUI::Cursor::Arrow);
    m_ibeamCursor = GUI::Cursor(GUI::Cursor::IBeam);
    m_waitCursor = GUI::Cursor(GUI::Cursor::Wait);
    m_crosshairCursor = GUI::Cursor(GUI::Cursor::Crosshair);
    m_waitArrowCursor = GUI::Cursor(GUI::Cursor::WaitArrow);
    m_sizeNWSECursor = GUI::Cursor(GUI::Cursor::SizeNWSE);
    m_sizeNESWCursor = GUI::Cursor(GUI::Cursor::SizeNESW);
    m_sizeWECursor = GUI::Cursor(GUI::Cursor::SizeWE);
    m_sizeAllCursor = GUI::Cursor(GUI::Cursor::SizeAll);
    m_noCursor = GUI::Cursor(GUI::Cursor::No);
    m_handCursor = GUI::Cursor(GUI::Cursor::Hand);
  }

  Plot::~Plot()
  {
    // Close all windows
    if (m_figures.size() > 0)
      m_figures.clear();

    // Cleanup the OpenGL context
    Magnum::GL::Context::makeCurrent(nullptr);
    m_context = std::nullopt;
    SDL_DestroyWindow(m_initWindow);
    SDL_GL_DeleteContext(m_glContext);

    // Finally quit SDL
    SDL_Quit();
  }

  bool Plot::setSwapInterval(SwapInterval interval)
  {
    if (SDL_GL_SetSwapInterval(interval) == -1)
    {
      std::cerr << "Failed to set swap interval" << std::endl;
      return false;
    }
    if (SDL_GL_GetSwapInterval() != interval)
    {
      std::cerr << "Swap interval ignored by driver" << std::endl;
      return false;
    }
    if (interval)
      m_isVsyncEnabled = true;
    else
      m_isVsyncEnabled = false;
    return true;
  }

  Plot::SwapInterval Plot::getSwapInterval() const
  {
    return static_cast<SwapInterval>(SDL_GL_GetSwapInterval());
  }

  void Plot::setMinimalLoopPeriod(unsigned int period)
  {
    m_minimalLoopPeriod = period;
  }

  unsigned int Plot::getMinimalLoopPeriod() const
  {
    return m_minimalLoopPeriod;
  }

  void Plot::show()
  {
    m_quit = false;

    // Draw all figures
    for (auto& [id, fig] : m_figures)
    {
      fig->makeCurrent();
      fig->clear();
      fig->drawContent();
      fig->swapBuffers();
    }

    // Show all figures
    for (auto& [id, fig] : m_figures)
      fig->setVisible(true);

    const Uint32 timeBefore =
      this->getMinimalLoopPeriod() ? SDL_GetTicks() : 0;

    while (!m_quit)
    {
      SDL_Event e;
      SDL_WaitEvent(&e);
      dispatchSDLEvent(std::move(e));

      SDL_Event additionalEvent;
      while (SDL_PollEvent(&additionalEvent))
        dispatchSDLEvent(std::move(additionalEvent));

      // If VSync is not enabled, delay to prevent CPU hogging (if set)
      if(getSwapInterval() == IMMEDIATE && getMinimalLoopPeriod())
      {
        const Uint32 loopTime = SDL_GetTicks() - timeBefore;
        if(loopTime < getMinimalLoopPeriod())
          SDL_Delay(getMinimalLoopPeriod() - loopTime);
      }

      // TODO: We can avoid some performance costs by tracking which figures
      // requiring redrawing
      for (auto& [id, fig] : m_figures)
      {
        if (fig->redraw())
        {
          fig->makeCurrent();
          fig->clear();
          fig->drawContent();
          fig->swapBuffers();
        }
      }
    } // End of main loop
  }

  void Plot::closeFigure(const FigureId& id)
  {
    m_figures.erase(m_figures.find(id));
  }

  void Plot::quit()
  {
    m_figures.clear();
    m_quit = true;
  }

  const SDL_GLContext& Plot::getGLContext() const
  {
    return m_glContext;
  }

  void Plot::dispatchSDLEvent(const SDL_Event& e)
  {
    switch (e.type)
    {
      case SDL_WINDOWEVENT:
      {
        switch (e.window.event)
        {
          case SDL_WINDOWEVENT_CLOSE:
          {
            closeFigure(e.window.windowID);
            break;
          }
        }
        break;
      }
      case SDL_MOUSEMOTION:
      if (m_figures.count(e.motion.windowID))
      {
        auto height = m_figures.at(e.motion.windowID)->getWindowSize().y();
        m_figures.at(e.motion.windowID)->handle(
            Backend::Event::MouseMotionEvent{
            e.motion.timestamp,
            e.motion.state,
            e.motion.x,
            height - e.motion.y,
            e.motion.xrel,
            -e.motion.yrel
            });
      }
      break;
      case SDL_MOUSEBUTTONUP:
      case SDL_MOUSEBUTTONDOWN:
      if (m_figures.count(e.button.windowID))
      {
        auto height = m_figures.at(e.motion.windowID)->getWindowSize().y();
        m_figures.at(e.button.windowID)->handle(
            Backend::Event::MouseButtonEvent{
            e.button.timestamp,
            Backend::Event::MouseEvent::Button(e.button.button),
            Backend::Event::MouseButtonEvent::ButtonState(e.button.state),
            Backend::Event::MouseButtonEvent::Clicks(e.button.clicks),
            e.button.x,
            height - e.button.y
            });
      }
      break;
      case SDL_MOUSEWHEEL:
      if (m_figures.count(e.wheel.windowID))
      {
        m_figures.at(e.wheel.windowID)->handle(
            Backend::Event::MouseWheelEvent{
            e.wheel.timestamp,
            e.wheel.x,
            e.wheel.y,
            Backend::Event::MouseWheelEvent::WheelDirection(e.wheel.direction)
            });
      }
      break;
      case SDL_QUIT:
      quit();
      return;
    }

    if (m_figures.size() == 0)
      quit();
  }

  Magnum::Math::Vector<Scalar>2<int> Plot::getMousePosition() const
  {
    int x, y;
    SDL_GetMouseState(&x, &y);
    return {x, y};
  }

  void Plot::setCursor(GUI::Cursor::SystemCursor c)
  {
    switch (c)
    {
      case GUI::Cursor::Arrow:
      SDL_SetCursor(m_arrowCursor.getHandle());
      break;
      case GUI::Cursor::IBeam:
      SDL_SetCursor(m_ibeamCursor.getHandle());
      break;
      case GUI::Cursor::Wait:
      SDL_SetCursor(m_waitCursor.getHandle());
      break;
      case GUI::Cursor::Crosshair:
      SDL_SetCursor(m_crosshairCursor.getHandle());
      break;
      case GUI::Cursor::WaitArrow:
      SDL_SetCursor(m_waitArrowCursor.getHandle());
      break;
      case GUI::Cursor::SizeNWSE:
      SDL_SetCursor(m_sizeNWSECursor.getHandle());
      break;
      case GUI::Cursor::SizeNESW:
      SDL_SetCursor(m_sizeNESWCursor.getHandle());
      break;
      case GUI::Cursor::SizeWE:
      SDL_SetCursor(m_sizeWECursor.getHandle());
      break;
      case GUI::Cursor::SizeAll:
      SDL_SetCursor(m_sizeAllCursor.getHandle());
      break;
      case GUI::Cursor::No:
      SDL_SetCursor(m_noCursor.getHandle());
      break;
      case GUI::Cursor::Hand:
      SDL_SetCursor(m_handCursor.getHandle());
      break;
    }
  }
}
