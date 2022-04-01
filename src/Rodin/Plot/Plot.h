/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_PLOT_H
#define RODIN_PLOT_PLOT_H

#include <map>
#include <memory>

#include <Magnum/Platform/GLContext.h>

#include "Rodin/Core/Common.h"

#include "Common.h"
#include "GUI/Cursor.h"
#include "Configuration.h"
#include "Artist/Figure.h"


namespace Rodin::Plot
{
  class Plot
  {
    private:
      using UPFigure = std::unique_ptr<Artist::Figure>;
      using FigureMap = std::map<FigureId, UPFigure>;

    public:
      /**
       * Represents the possible values of the swap interval for objects of the
       * Artist::Figure class.
       * @see setSwapInterval
       */
      enum SwapInterval
      {
        /**
         * Updates synchronized with the vertical retrace, except that if the
         * vertical retrace for the current frame was missed the buffers are
         * swapped immediately.
         */
        ADAPTIVE = -1,

        /**
         * Immediate updates.
         */
        IMMEDIATE = 0,

        /**
         * Updates synchronized with the vertical retrace.
         */
        VSYNC = 1
      };

      Plot(const Plot&)   = delete;
      Plot(Plot&&)        = delete;
      void operator=(const Plot&)   = delete;
      void operator=(Plot&&)        = delete;

      Plot();
      Plot(const Configuration& configuration);
      ~Plot();

      /**
       * Creates a new figure managed by the Plot object.
       * @param[in] args Arguments to be passed to the Figure constructor
       * @returns Reference to the Artist::Figure instance.
       * @see Artist::Figure
       */
      template <class ... Args>
      Artist::Figure& figure(Args&&... args)
      {
        auto fig = std::make_unique<Artist::Figure>(
            *this, std::forward<Args>(args)...);
        auto id = fig->getId();
        m_figures.insert({id, std::move(fig)});
        return *m_figures.at(id);
      }

      /**
       * Creates a new figure managed by the Plot object.
       * @param[in] args Arguments to be passed to the Figure constructor
       * @returns Const reference to the Artist::Figure instance.
       * @see Artist::Figure
       */
      template <class ... Args>
      const Artist::Figure& figure(Args&&... args) const
      {
        auto fig = std::make_unique<Artist::Figure>(
            *this, std::forward<Args>(args)...);
        auto id = fig->getId();
        m_figures.insert({id, std::move(fig)});
        return *m_figures.at(id);
      }

      /**
       * Sets the global swap interval for each figure managed by the Plot
       * object.
       *
       * @param[in] interval Swap interval
       * @returns `true` if the operation was successful, `false` otherwise.
       * @note Default value is drived dependent.
       */
      bool setSwapInterval(SwapInterval interval);

      /**
       * Queries the current value of the swap interval.
       * @returns Swap interval
       */
      SwapInterval getSwapInterval() const;

      void setMinimalLoopPeriod(unsigned int period);

      unsigned int getMinimalLoopPeriod() const;

      /**
       * Displays all the figures and enters the event handling loop, which
       * blocks execution.
       */
      void show();

      /**
       * @returns The number of figures managed by the the Plot object.
       */
      size_t count() const;


      /**
       * @returns The underlying OpenGL context.
       */
      const SDL_GLContext& getGLContext() const;

      Eigen::Array2<int> getMousePosition() const;

      /**
       * Sets the cursor 
       */
      void setCursor(GUI::Cursor::SystemCursor c);

    private:
      void dispatchSDLEvent(const SDL_Event& e);
      void closeFigure(const FigureId& id);
      void quit();

      int     m_argc;
      char**  m_argv;

      unsigned int    m_minimalLoopPeriod;
      bool            m_isVsyncEnabled,
                      m_quit;

      FigureMap       m_figures;
      WindowHandle                                  m_initWindow;
      SDL_GLContext                                 m_glContext;
      std::optional<Magnum::Platform::GLContext>    m_context;

      GUI::Cursor  m_arrowCursor,
                   m_ibeamCursor,
                   m_waitCursor,
                   m_crosshairCursor,
                   m_waitArrowCursor,
                   m_sizeNWSECursor,
                   m_sizeNESWCursor,
                   m_sizeWECursor,
                   m_sizeAllCursor,
                   m_noCursor,
                   m_handCursor;
  };
}

#endif
