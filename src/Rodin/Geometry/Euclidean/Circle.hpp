#ifndef RODIN_GEOMETRY_EUCLIDEAN_CIRCLE_IPP
#define RODIN_GEOMETRY_EUCLIDEAN_CIRCLE_IPP

#include <cmath>

#include "Point2D.h"
#include "Line2D.h"
#include "LineSegment2D.h"

#include "Circle.h"

namespace Rodin::Geometry::Euclidean
{
  template <class T>
  constexpr
  Circle<T>::Circle(const Point2D<T>& center, const T& radius)
    : m_center(center), m_radius(radius)
  {
    assert(radius > T(0));
  }

  template <class T>
  inline
  constexpr
  T Circle<T>::operator()(const Point2D<T>& p) const
  {
    return (p.x() - center().x()) * (p.x() - center().x())
      + (p.y() - center().y()) * (p.y() - center().y()) - radius() * radius();
  }

  template <class T>
  inline
  constexpr
  Point2D<T> Circle<T>::operator()(const Math::Rad& angle) const
  {
    return
      m_center + m_radius * Point2D<T>(
          std::cos(static_cast<Real>(angle)), std::sin(static_cast<Real>(angle)));
  }

  template <class T>
  inline
  constexpr
  T Circle<T>::radius() const
  {
    return m_radius;
  }

  template <class T>
  inline
  constexpr
  Circle<T>& Circle<T>::setRadius(T radius)
  {
    m_radius = radius;
    return *this;
  }

  template <class T>
  inline
  constexpr
  Point2D<T> Circle<T>::center() const
  {
    return m_center;
  }

  template <class T>
  inline
  constexpr
  Circle<T>& Circle<T>::setCenter(const Point2D<T>& center)
  {
    m_center = center;
    return *this;
  }

  template <class T>
  inline
  constexpr
  Line2D<T> Circle<T>::tangent(const Math::Rad& angle)
  {
    auto p = operator()(angle);
    return Line2D<T>(
        p.x() - center().x(),
        p.y() - center().y(),
        p.x() * p.x() + p.y() * p.y() - center().x() * p.x() - center().y() * p.y()
        );
  }
}

#endif
