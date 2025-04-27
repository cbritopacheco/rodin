/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <Magnum/EigenIntegration/Integration.h>

#include "Rodin/Plot/Artist/Figure.h"

#include "BaseAxes.h"

namespace Rodin::Plot::Backend::Bases
{
  BaseAxes::BaseAxes(
    Artist::Figure& figure,
    const Magnum::Math::Vector2<Integer>& bottomLeft,
    const Magnum::Math::Vector2<Integer>& size,
    Boolean frameEnabled)
    : m_figure(figure),
      m_bottomLeft(bottomLeft),
      m_size(size),
      m_frameEnabled(frameEnabled)
  {}

  const Magnum::Math::Vector2<Integer>& BaseAxes::getBottomLeft() const
  {
    return m_bottomLeft;
  }

  const Magnum::Math::Vector2<Integer>& BaseAxes::getSize() const
  {
    return m_size;
  }

  Artist::Figure& BaseAxes::getFigure()
  {
    return m_figure;
  }

  const Artist::Figure& BaseAxes::getFigure() const
  {
    return m_figure;
  }

  Boolean BaseAxes::isFrameEnabled() const
  {
    return m_frameEnabled;
  }

  BaseAxes::XLimits BaseAxes::getXLimits() const
  {
    return m_xlim;
  }

  BaseAxes::YLimits BaseAxes::getYLimits() const
  {
   return m_ylim;
  }

  BaseAxes& BaseAxes::enableFrame(Boolean v)
  {
   m_frameEnabled = v;
   return *this;
  }

  BaseAxes& BaseAxes::setXLimits(const XLimits& xlim)
  {
    assert(xlim.left < xlim.right);
    m_xlim = xlim;
    return *this;
  }

  BaseAxes& BaseAxes::setYLimits(const YLimits& ylim)
  {
    assert(ylim.bottom < ylim.top);
    m_ylim = ylim;
    return *this;
  }

  Geometry::Euclidean::Rectangle<Integer> BaseAxes::getBoundingBox() const
  {
    return Geometry::Euclidean::Rectangle<Integer>(
        { getBottomLeft().x(), getBottomLeft().y() },
        { getBottomLeft().x() + getSize().x(), getBottomLeft().y() + getSize().y() });
  }
}
