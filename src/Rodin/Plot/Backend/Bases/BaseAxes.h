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

#include "Rodin/Types.h"
#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Event.h"
#include "Rodin/Geometry/Euclidean/Rectangle.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseAxes
  {
    public:
      struct XLimits
      {
        Float left, right;
      };

      struct YLimits
      {
        Float bottom, top;
      };

      BaseAxes(
          Artist::Figure& fig,
          const Magnum::Math::Vector<Real>2<Integer>& bottomLeft,
          const Magnum::Math::Vector<Real>2<Integer>& size,
          Boolean frameEnabled);

      virtual ~BaseAxes() = default;

      virtual void drawContent() = 0;

      const Magnum::Math::Vector<Real>2<Integer>& getSize() const;

      const Magnum::Math::Vector<Real>2<Integer>& getBottomLeft() const;

      Artist::Figure& getFigure();
      const Artist::Figure& getFigure() const;

      bool isFrameEnabled() const;

      BaseAxes& enableFrame(Boolean v);

      XLimits getXLimits() const;

      YLimits getYLimits() const;

      BaseAxes& setXLimits(const XLimits& xlim);
      BaseAxes& setYLimits(const YLimits& ylim);

      Geometry::Euclidean::Rectangle<Integer> getBoundingBox() const;

      virtual void handle(const Backend::Event::MouseMotionEvent& e) = 0;
      virtual void handle(const Backend::Event::MouseButtonEvent& e) = 0;
      virtual void handle(const Backend::Event::MouseWheelEvent& e) = 0;

    private:
      Artist::Figure&     m_figure;

      Magnum::Math::Vector<Real>2<Integer> m_bottomLeft;
      Magnum::Math::Vector<Real>2<Integer> m_size;

      Boolean     m_frameEnabled;
      XLimits     m_xlim;
      YLimits     m_ylim;
  };
}

#endif
