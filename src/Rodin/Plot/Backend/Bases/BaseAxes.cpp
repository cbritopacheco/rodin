/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>

#include "Rodin/Plot/Artist/Figure.h"

#include "BaseAxes.h"

namespace Rodin::Plot::Backend::Bases
{
  BaseAxes::BaseAxes(
      Artist::Figure& figure,
      Eigen::Array2<int> bottomLeft,
      Eigen::Array2<int> size,
      bool frameEnabled)
    : m_figure(figure),
      m_bottomLeft(bottomLeft),
      m_size(size),
      m_frameEnabled(frameEnabled)
  {}

  Eigen::Array2<int> BaseAxes::getBottomLeft() const
  {
    return m_bottomLeft;
  }

  Eigen::Array2<int> BaseAxes::getSize() const
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

  bool BaseAxes::isFrameEnabled() const
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

  BaseAxes& BaseAxes::enableFrame(bool v)
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

  Core::Geometry::Rectangle<int> BaseAxes::getBoundingBox() const
  {
    return Core::Geometry::Rectangle<int>(
        getBottomLeft().matrix(),
        {
          getBottomLeft().x() + getSize().x(),
          getBottomLeft().y() + getSize().y()
        });
  }
}
