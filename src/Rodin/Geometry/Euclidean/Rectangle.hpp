#ifndef RODIN_GEOMETRY_EUCLIDEAN_RECTANGLE_IPP
#define RODIN_GEOMETRY_EUCLIDEAN_RECTANGLE_IPP

#include "Point2D.h"
#include "Rectangle.h"

namespace Rodin::Geometry::Euclidean
{
  template <class T>
  constexpr
  Rectangle<T>::Rectangle(const Point2D<T>& bottomLeft, const Point2D<T>& topRight)
    : m_bottomLeft(bottomLeft), m_topRight(topRight)
  {}

  template <class T>
  inline
  constexpr
  bool Rectangle<T>::contains(const Point2D<T>& p) const
  {
    return
      (getCorner<BottomLeft>().x() <= p.x()) &&
      (getCorner<BottomLeft>().y() <= p.y()) &&
      (p.x() <= getCorner<TopRight>().x()) &&
      (p.y() <= getCorner<TopRight>().y());
  }

  template <class T>
  template <typename Rectangle<T>::Corner c>
  inline
  constexpr
  Point2D<T> Rectangle<T>::getCorner() const
  {
    if constexpr (c == BottomLeft)
    {
      return m_bottomLeft;
    }
    else if constexpr (c == BottomRight)
    {
      return m_bottomLeft + Eigen::Vector2<T>(m_topRight.x() - m_bottomLeft.x(), 0);
    }
    else if constexpr (c == TopLeft)
    {
      return m_bottomLeft + Eigen::Vector2<T>(0, m_topRight.y() - m_bottomLeft.y());
    }
    else if constexpr (c == TopRight)
    {
      return m_topRight;
    }
    else
    {
      // Well if we arrived here then something is horribly wrong
      assert(false);
    }
  }

  template <class T>
  inline
  constexpr
  T Rectangle<T>::height() const
  {
    return m_topRight.y() - m_bottomLeft.y();
  }

  template <class T>
  inline
  constexpr
  T Rectangle<T>::width() const
  {
    return m_topRight.x() - m_bottomLeft.x();
  }

  template <class T>
  inline
  constexpr
  T Rectangle<T>::area() const
  {
    return width() * height();
  }
}

#endif
