/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LineSegment.hpp
 * @brief Line segment operations in a compile-time Euclidean dimension.
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_LINESEGMENT_HPP
#define RODIN_GEOMETRY_EUCLIDEAN_LINESEGMENT_HPP

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "Point.h"

#include "LineSegment.h"
#include "Line.h"

namespace Rodin::Geometry::Euclidean
{

  template <size_t Dimension, class T>
  constexpr LineSegment<Dimension, T>::LineSegment(
    const Point<Dimension, T>& start, const Point<Dimension, T>& end)
    : m_start(start),
      m_end(end)
  {
    if (start == end)
      Alert::MemberFunctionException(*this, __func__)
        << "LineSegment requires distinct endpoints" << Alert::Raise;
  }

  template <size_t Dimension, class T>
  inline constexpr Point<Dimension, T> LineSegment<Dimension, T>::operator()(
    const T& t) const
  {
    return start() + t * (end() - start());
  }

  template <size_t Dimension, class T>
  inline constexpr T LineSegment<Dimension, T>::length() const
  {
    if constexpr (Dimension == 0)
      return T(0);
    else if constexpr (Dimension == 1)
      return std::abs(m_end[0] - m_start[0]);
    else
      return (m_end - m_start).norm();
  }

  template <size_t Dimension, class T>
  inline constexpr typename Point<Dimension, T>::Vector
  LineSegment<Dimension, T>::direction() const
  {
    if constexpr (Dimension == 0)
      return Point<Dimension, T>::Vector::Zero();
    else if constexpr (Dimension == 1)
    {
      typename Point<Dimension, T>::Vector result;
      result[0] = m_end[0] > m_start[0] ? T(1) : T(-1);
      return result;
    }
    else
      return (m_end - m_start).normalized();
  }

  template <size_t Dimension, class T>
  inline constexpr Point<Dimension, T> LineSegment<Dimension, T>::start() const
  {
    return m_start;
  }

  template <size_t Dimension, class T>
  inline constexpr Point<Dimension, T> LineSegment<Dimension, T>::end() const
  {
    return m_end;
  }

  template <size_t Dimension, class T>
  inline constexpr LineSegment<Dimension, T>& LineSegment<Dimension, T>::setStart(
    const Point<Dimension, T>& start)
  {
    if (start == m_end)
      Alert::MemberFunctionException(*this, __func__)
        << "LineSegment requires distinct endpoints" << Alert::Raise;
    m_start = start;
    return *this;
  }

  template <size_t Dimension, class T>
  inline constexpr LineSegment<Dimension, T>& LineSegment<Dimension, T>::setEnd(
    const Point<Dimension, T>& end)
  {
    if (end == m_start)
      Alert::MemberFunctionException(*this, __func__)
        << "LineSegment requires distinct endpoints" << Alert::Raise;
    m_end = end;
    return *this;
  }

  template <size_t Dimension, class T>
  inline constexpr LineSegment<Dimension, T> LineSegment<Dimension, T>::reverse() const
  {
    return LineSegment<Dimension, T>(end(), start());
  }

  template <size_t Dimension, class T>
  inline constexpr Optional<Point<Dimension, T>> LineSegment<Dimension, T>::intersect(
    const LineSegment<Dimension, T>& other) const
  {
    const Line<Dimension, T> supporting(*this), otherSupporting(other);
    const auto point = supporting.intersect(otherSupporting);
    if (point)
    {
      const auto u = m_end - m_start, v = other.m_end - other.m_start;
      const T s = u.dot(*point - m_start) / u.squaredNorm();
      const T t = v.dot(*point - other.m_start) / v.squaredNorm();
      if (T(0) <= s && s <= T(1) && T(0) <= t && t <= T(1))
        return point;
      return {};
    }
    // Parallel segments can meet at a unique shared endpoint. A positive-length
    // overlap has no unique intersection point and returns no point.
    const auto u = m_end - m_start;
    const T scale = std::max(
      {T(1), m_start.norm(), m_end.norm(), other.m_start.norm(), other.m_end.norm()});
    const T tolerance = T(Dimension) * std::numeric_limits<T>::epsilon() * scale;
    if (supporting.distance(other.m_start) > tolerance ||
      supporting.distance(other.m_end) > tolerance)
      return {};
    const T a = u.dot(other.m_start - m_start) / u.squaredNorm();
    const T b = u.dot(other.m_end - m_start) / u.squaredNorm();
    const T lower = std::max(T(0), std::min(a, b));
    const T upper = std::min(T(1), std::max(a, b));
    if (lower == upper)
      return operator()(lower);
    return {};
  }
}

#endif
