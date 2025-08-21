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

#include "Common.h"
#include "GUI/Cursor.h"
#include "Configuration.h"
#include "Artist/Figure.h"


namespace Rodin::Plot
{
  /**
   * @defgroup RodinPlot Plotting and Visualization Module
   * @brief Interactive plotting and visualization capabilities for finite element data.
   *
   * The Plot module provides comprehensive visualization tools for finite element
   * analysis results, including mesh visualization, function plotting, and interactive
   * exploration of numerical solutions. Built on top of the Magnum graphics library,
   * it offers high-performance rendering suitable for large-scale simulations.
   *
   * ## Key Features
   * - **Mesh Visualization**: Wireframe and solid mesh rendering
   * - **Function Plotting**: Scalar and vector field visualization
   * - **Interactive Controls**: Mouse and keyboard interaction for exploration
   * - **Multiple Backends**: Support for different rendering backends
   * - **Export Capabilities**: Save visualizations to various image formats
   */

  /**
   * @ingroup RodinPlot
   * @brief Main plotting interface for finite element visualization.
   *
   * The Plot class serves as the central interface for creating and managing
   * visualizations of finite element data. It provides a high-level API for
   * rendering meshes, functions, and analysis results with interactive controls.
   *
   * ## Usage Example
   * ```cpp
   * Plot plot;
   * auto figure = plot.figure();
   * figure.plot(mesh);
   * figure.plot(solution);
   * plot.show();
   * ```
   *
   * ## Technical Details
   * - Uses OpenGL for hardware-accelerated rendering
   * - Supports multiple concurrent figure windows
   * - Provides customizable rendering options and visual styles
   * - Integrates with Rodin's finite element data structures
   */
  class Plot
  {
    private:
      using UPFigure = std::unique_ptr<Artist::Figure>;
      using FigureMap = std::map<FigureId, UPFigure>;

    public:
      /**
       * @brief Buffer swap interval modes for controlling rendering synchronization.
       *
       * This enumeration defines the synchronization behavior between buffer
       * swaps and the display's vertical refresh rate, allowing control over
       * frame rate and visual quality trade-offs.
       */
      enum SwapInterval
      {
        /**
         * @brief Adaptive vertical synchronization.
         *
         * Updates synchronized with the vertical retrace, except that if the
         * vertical retrace for the current frame was missed the buffers are
         * swapped immediately. This provides smooth rendering while avoiding
         * input lag.
         */
        ADAPTIVE = -1,

        /**
         * @brief Immediate buffer swaps.
         *
         * Immediate updates without synchronization. This provides the highest
         * frame rate but may result in screen tearing.
         */
        IMMEDIATE = 0,

        /**
         * @brief Vertical synchronization.
         *
         * Updates synchronized with the vertical retrace. This eliminates
         * screen tearing but may introduce input lag and limits frame rate
         * to the display refresh rate.
         */
        VSYNC = 1
      };

      Plot(const Plot&)  = delete;

      Plot(Plot&&)      = delete;

      void operator=(const Plot&)  = delete;

      void operator=(Plot&&)      = delete;

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
        auto fig = std::make_unique<Artist::Figure>(*this, std::forward<Args>(args)...);
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
        auto fig = std::make_unique<Artist::Figure>(*this, std::forward<Args>(args)...);
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

      Magnum::Math::Vector2<int> getMousePosition() const;

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

      unsigned int  m_minimalLoopPeriod;
      bool          m_isVsyncEnabled,
                    m_quit;

      FigureMap       m_figures;
      WindowHandle    m_initWindow;
      SDL_GLContext   m_glContext;
      Optional<Magnum::Platform::GLContext>  m_context;

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
