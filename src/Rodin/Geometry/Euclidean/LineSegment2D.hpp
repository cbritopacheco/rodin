#ifndef RODIN_GEOMETRY_EUCLIDEAN_SEGMENT2D_IPP
#define RODIN_GEOMETRY_EUCLIDEAN_SEGMENT2D_IPP

#include <Eigen/Core>

#include "Point2D.h"
#include "Circle.h"

#include "LineSegment2D.h"

namespace Rodin::Geometry::Euclidean
{

  template <class T>
  constexpr
  LineSegment2D<T>::LineSegment2D(
      const Point2D<T>& start, const Point2D<T>& end)
    : m_start(start), m_end(end)
  {
    assert(start != end);
  }

  template <class T>
  inline
  constexpr
  Point2D<T> LineSegment2D<T>::operator()(const T& t) const
  {
    return start() + t * (end() - start());
  }

  template <class T>
  inline
  constexpr
  T LineSegment2D<T>::length() const
  {
    return (m_end - m_start).norm();
  }

  template <class T>
  inline
  constexpr
  Eigen::Vector2<T> LineSegment2D<T>::direction() const
  {
    return (m_end - m_start).normalized();
  }

  template <class T>
  inline
  constexpr
  Point2D<T> LineSegment2D<T>::start() const
  {
    return m_start;
  }

  template <class T>
  inline
  constexpr
  Point2D<T> LineSegment2D<T>::end() const
  {
    return m_end;
  }

  template <class T>
  inline
  constexpr
  LineSegment2D<T>& LineSegment2D<T>::setStart(const Point2D<T>& start)
  {
    m_start = start;
    return *this;
  }

  template <class T>
  inline
  constexpr
  LineSegment2D<T>& LineSegment2D<T>::setEnd(const Point2D<T>& end)
  {
    m_end = end;
    return *this;
  }

  template <class T>
  inline
  constexpr
  LineSegment2D<T> LineSegment2D<T>::reverse() const
  {
    return LineSegment2D<T>(end(), start());
  }

  template <class T>
  inline
  constexpr
  std::optional<Point2D<T>> LineSegment2D<T>::intersect(const
      LineSegment2D<T>& other)
  const
  {
    auto det = (end().x() - start().x()) * (other.start().y() - other.end().y()) -
      (end().y() - start().y()) * (other.start().x() - other.end().x());

    if (det == T{0})
      return {};
    else
    {
      auto s0 = (other.start().y() - other.end().y()) * (other.start().x() - start().x())
            + (other.end().x() - other.start().x()) * (other.start().y() - start().y());
      s0 /= det;

      auto t0 = (start().y() - end().y()) * (other.start().x() - start().x())
            + (end().x() - start().x()) * (other.start().y() - start().y());
      t0 /= det;

      if (T{0} <= s0 && s0 <= T{1} && T{0} <= t0 && t0 <= T{1})
        return other(t0);
      else
        return {};
    }
  }
}

#endif
