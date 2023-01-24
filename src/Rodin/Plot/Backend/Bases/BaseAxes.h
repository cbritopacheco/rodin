/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_BASES_BASEAXES_H
#define RODIN_PLOT_BACKEND_BASES_BASEAXES_H

#include <memory>
#include <vector>

#include <Magnum/Array.h>

#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Event.h"
#include "Rodin/Plot/Geometry/Rectangle.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseAxes
  {
   public:
    struct XLimits
    {
      float left, right;
    };

    struct YLimits
    {
      float bottom, top;
    };

    BaseAxes(
       Artist::Figure& fig,
       Magnum::Math::Vector2<int> bottomLeft,
       Magnum::Math::Vector2<int> size,
       bool frameEnabled);

    virtual ~BaseAxes() = default;

    virtual void drawContent() = 0;

    Magnum::Math::Vector2<int> getSize() const;

    Magnum::Math::Vector2<int> getBottomLeft() const;

    Artist::Figure& getFigure();
    const Artist::Figure& getFigure() const;

    bool isFrameEnabled() const;

    BaseAxes& enableFrame(bool v);

    XLimits getXLimits() const;

    YLimits getYLimits() const;

    BaseAxes& setXLimits(const XLimits& xlim);
    BaseAxes& setYLimits(const YLimits& ylim);

    Geometry::Rectangle<int> getBoundingBox() const;

    virtual void handle(const Backend::Event::MouseMotionEvent& e) = 0;
    virtual void handle(const Backend::Event::MouseButtonEvent& e) = 0;
    virtual void handle(const Backend::Event::MouseWheelEvent& e) = 0;
   private:
    Artist::Figure&     m_figure;

    Magnum::Math::Vector2<int> m_bottomLeft;
    Magnum::Math::Vector2<int> m_size;

    bool       m_frameEnabled;
    XLimits     m_xlim;
    YLimits     m_ylim;
  };
}

#endif
