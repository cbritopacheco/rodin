#ifndef RODIN_CORE_GEOMETRY_POINT2D_IPP
#define RODIN_CORE_GEOMETRY_POINT2D_IPP

#include "Line2D.h"
#include "Circle.h"
#include "LineSegment2D.h"

#include "Point2D.h"

namespace Rodin::Plot::Geometry
{
  template <class T>
  inline
  constexpr
  std::optional<LineSegment2D<T>> Point2D<T>::connect(const Point2D<T>& other) const
  {
    if (*this == other)
      return {};
    else
      return LineSegment2D<T>(*this, other);
  }

  template <class T>
  inline
  constexpr
  std::optional<LineSegment2D<T>> Point2D<T>::connect(const Line2D<T>& line) const
  {
    assert(line.a() != T{0} || line.b() != T{0});

    if (line(*this) > T{0} || line(*this) < T{0}) // Point is not touching the line
    {
      return LineSegment2D<T>(
          *this,
          Point2D<T>({
            line.b() * (line.b() * Magnum::Math::Vector2<T>::x() - line.a() *
                Magnum::Math::Vector2<T>::y()) + line.a() * line.c(),
            line.a() * (-line.b() * Magnum::Math::Vector2<T>::x() + line.a() *
                Magnum::Math::Vector2<T>::y()) + line.b() * line.c()
            }) / (line.a() * line.a() + line.b() * line.b()));
    }
    else
    {
      return {};
    }
  }

  template <class T>
  inline
  constexpr
  std::optional<LineSegment2D<T>> Point2D<T>::connect(const Circle<T>& other) const
  {
    auto v = other(*this);
    if (v < T{0})
    {
      auto centerToPoint = other.center().connect(*this);
      return LineSegment2D<T>(
          *this,
          (other.radius() - centerToPoint.length()) * centerToPoint.direction());
    }
    else if (v > T{0})
    {
      auto centerToPoint = other.center().connect(*this);
      return LineSegment2D<T>(
          *this,
          other.center() + other.radius() * centerToPoint.direction());
    }
    else
    {
      return {};
    }
  }
}

#endif
