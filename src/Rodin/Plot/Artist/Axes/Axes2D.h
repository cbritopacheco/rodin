/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_ARTIST_AXES_AXES2D_H
#define RODIN_PLOT_ARTIST_AXES_AXES2D_H

#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Bases/BaseAxes.h"
#include "Rodin/Plot/Backend/Bases/BaseTopLevelArtist2D.h"
#include "Rodin/Plot/Backend/Renderer/Camera2D.h"

namespace Rodin::Plot::Artist::Axes
{
  class Axes2D : public Backend::Bases::BaseAxes, public Backend::Bases::BaseTopLevelArtist2D
  {
    public:
      /**
      * Construct Axes2D object.
      * @param figRef Figure reference.
      * @param bottomLeft Bottom left coordinate.
      * @param size Width and height dimensions
      * @param frameEnable Enable the drawing of the axes frames
      */
      Axes2D(
         Figure& figRef,
         const Magnum::Math::Vector2<Integer>& bottomLeft,
         const Magnum::Math::Vector2<Integer>& size,
         bool frameEnabled = true);

      virtual ~Axes2D() = default;

      /**
      * Plots `y` versus `x` as lines and/or markers.
      * @param x Horizontal coordinates of data points
      * @param y Vertical coordinates of data points
      * @param Other arguments to be passed to the Artist::Line2D constructor
      * @returns Line artist for the manipulation of the line
      */
      template <class T, class ... Args>
      Lines::Line2D& plot(const Eigen::ArrayX<T>& x, const Eigen::ArrayX<T>& y, Args&&... args);

      virtual void drawContent() override;
      virtual void handle(const Backend::Event::MouseMotionEvent& e) override;
      virtual void handle(const Backend::Event::MouseButtonEvent& e) override;
      virtual void handle(const Backend::Event::MouseWheelEvent& e) override;

      virtual Backend::Renderer::DrawableGroup2D& getDrawableGroup() override;
      virtual const Backend::Renderer::DrawableGroup2D& getDrawableGroup() const override;

      /**
      * @returns The 2D object handle for manipulating the object coordinates
      * in 2D space.
      */
      virtual Backend::Renderer::Object2D& getObject2D() override;

   private:
    // Renderer objects
    Backend::Renderer::Scene2D   m_scene;
    Backend::Renderer::Camera2D  m_camera;

    Backend::Renderer::DrawableGroup2D    m_drawables;

    // XLimits m_xlim;
    // YLimits m_ylim;
  };
}

#include "Axes2D.hpp"

#endif
