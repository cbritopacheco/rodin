/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_ARTIST_FIGURE_H
#define RODIN_PLOT_ARTIST_FIGURE_H

#include <optional>
#include <functional>

#include "Rodin/Plot/Common.h"
#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Renderer/Common.h"
#include "Rodin/Plot/Backend/Renderer/Camera2D.h"
#include "Rodin/Plot/Backend/Renderer/Drawables/Frame.h"
#include "Rodin/Plot/Backend/Bases/BaseTopLevelArtist2D.h"

#include "Rodin/Plot/Backend/Bases/BaseFigure.h"

namespace Rodin::Plot::Artist
{
  class Figure : public Backend::Bases::BaseFigure, public Backend::Bases::BaseTopLevelArtist2D
  {
    struct AxesWrapper
    {
      std::unique_ptr<Backend::Bases::BaseAxes> axes;
      std::reference_wrapper<Backend::Renderer::Drawables::Frame> frame;
    };

    public:
      using Backend::Bases::BaseFigure::handle;

      /**
       * Constructs a Figure object with a window size of 800x600 pixels.
       * @param pltRef Plot reference
       */
      Figure(Plot& pltRef);

      /**
       * Constructs a Figure object.
       * @param pltRef Plot reference
       * @param width Width in pixels of the window
       * @param Height Height in pixels of the window
       */
      Figure(Plot& pltRef, int width, int height);

      /**
       * Constructs a Figure object.
       * @param pltRef Plot reference.
       * @param width Width in pixels of the window.
       * @param Height Height in pixels of the window.
       * @param title Title of the window.
       */
      Figure(Plot& pltRef, int width, int height, const std::string& title);

      ~Figure() = default;

      /**
       * Add an axes to the figure.
       * @tparam AxesT Type of the Axes which will be added. Default is
       *  Rodin::Plot::Axes2D.
       * @tparam Args Arguments to forward to the AxesT constructor.
       * @returns A reference to the added Axes.
       * @warning The returned reference will be invalidated once the figure
       * object is destroyed.
       */
      template <class AxesT = Axes::Axes2D, class ... Args>
        AxesT& addAxes(Args&&... args);

      /**
       * Add an axes to the figure.
       * @tparam AxesT Type of the Axes which will be added. Default is
       *  Rodin::Plot::Axes2D.
       * @warning The returned reference will be invalidated once the figure
       * object is destroyed.
       */
      template <class AxesT = Axes::Axes2D>
        AxesT& addAxes();

      /**
       * Draws the figure's axes to its frame buffer.
       */
      void drawContent() override;

      /**
       * Clears the figure's frame buffer.
       */
      void clear();

      void handle(const Backend::Event::MouseMotionEvent& e) override;
      void handle(const Backend::Event::MouseButtonEvent& e) override;
      void handle(const Backend::Event::MouseWheelEvent& e) override;

      Plot& getPlot();

      /**
       * Returns the owning Plot reference.
       */
      const Plot& getPlot() const;

      virtual Backend::Renderer::DrawableGroup2D& getDrawableGroup() override;
      virtual const Backend::Renderer::DrawableGroup2D& getDrawableGroup() const override;

      /**
       * Returns the 2D object handle.
       */
      Backend::Renderer::Object2D& getObject2D() override;

      /**
       * Indicates if the figure's contents require re-drawing.
       * @retval true If requires redrawing
       * @retval false If does not require redraw
       */
      bool redraw() const;

    private:
      static unsigned int s_figureCount;

      // General references
      Plot&                          m_pltRef;

      // Renderer objects
      Backend::Renderer::Scene2D            m_scene;
      Backend::Renderer::Camera2D           m_camera;
      Backend::Renderer::DrawableGroup2D      m_drawables;

      // Managed objects
      std::vector<AxesWrapper>       m_axes;

      bool m_redraw;
  };
}

#include "Figure.hpp"

#endif
