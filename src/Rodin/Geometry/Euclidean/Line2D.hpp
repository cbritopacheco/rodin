#ifndef RODIN_GEOMETRY_EUCLIDEAN_LINE2D_IPP
#define RODIN_GEOMETRY_EUCLIDEAN_LINE2D_IPP

#include "Point2D.h"
#include "Circle.h"
#include "LineSegment2D.h"

#include "Line2D.h"

namespace Rodin::Geometry::Euclidean
{
  template <class T>
  constexpr
  Line2D<T>::Line2D(T a, T b, T c)
    : m_a(a), m_b(b), m_c(c)
  {
    assert(m_a != T{0} || m_b != T{0});
  }

  template <class T>
  constexpr
  Line2D<T>::Line2D(const Point2D<T>& p, const Point2D<T>& q)
    : m_a(p.y() - q.y()),
      m_b(q.x() - p.x()),
      m_c(q.x() * p.y() - p.x() * q.y())
  {
    assert(p != q);
    assert(m_a != T{0} || m_b != T{0});
  }

  template <class T>
  constexpr
  Line2D<T>::Line2D(std::initializer_list<T> list)
    : Line2D<T>(*list.begin(), *(list.begin() + 1), *(list.begin() + 2))
  {
    assert(list.size() == 3);
  }

  template <class T>
  constexpr
  Line2D<T>::Line2D(const LineSegment2D<T>& ls)
    : Line2D<T>(ls.start(), ls.end())
  {}

  template <class T>
  inline
  constexpr
  std::optional<T> Line2D<T>::operator()(T x) const
  {
    if (b() == T{0})
      return {};
    else
      return (c() - a() * x) / b();
  };

  template <class T>
  inline
  constexpr
  T Line2D<T>::operator()(const Point2D<T>& p) const
  {
    return a() * p.x() + b() * p.y() - c();
  }

  template <class T>
  inline
  constexpr
  T Line2D<T>::a() const
  {
    return m_a;
  }

  template <class T>
  inline
  constexpr
  T Line2D<T>::b() const
  {
    return m_b;
  }

  template <class T>
  constexpr
  T Line2D<T>::c() const
  {
    return m_c;
  }

  template <class T>
  inline
  constexpr
  std::optional<T> Line2D<T>::slope() const
  {
    if (b() == T{0})
      return {};
    else
      return -a() / b();
  }

  template <class T>
  inline
  constexpr
  std::optional<T> Line2D<T>::xIntercept() const
  {
    if (a() == T{0})
      return {};
    else
      return c() / a();
  }

  template <class T>
  inline
  constexpr
  std::optional<T> Line2D<T>::yIntercept() const
  {
    if (b() == T{0})
      return {};
    else
      return c() / b();
  }

  template <class T>
  inline
  constexpr
  bool Line2D<T>::isVertical() const
  {
    return b() == T{0};
  }

  template <class T>
  inline
  constexpr
  bool Line2D<T>::isHorizontal() const
  {
    return a() == T{0};
  }

  template <class T>
  inline
  constexpr
  std::optional<Point2D<T>> Line2D<T>::intersect(const Line2D<T>& other) const
  {
    auto det = a() * other.b() - other.a() * b();
    if (det == T{0})
    {
      return {};
    } else
    {
      return Point2D<T>({
          c() * other.b() - b() * other.c(),
          a() * other.c() - other.a() * c()}) / det;
    }
  }

  template <class T>
  inline
  constexpr
  T Line2D<T>::distance(const Line2D<T>& other) const
  {
    if (isVertical() && other.isVertical())
      return Rodin::Math::abs<T>(
          *xIntercept() - *other.xIntercept());
    else if (!isVertical() && !other.isVertical())
      return Rodin::Math::abs<T>(
          *yIntercept() - *other.yIntercept()) / Rodin::Math::sqrt(*slope() *
          *slope() + T{1});
    else
      return T{0};
  }

  template <class T>
  inline
  constexpr
  std::variant<std::nullopt_t, Point2D<T>, std::pair<Point2D<T>, Point2D<T>>>
  Line2D<T>::intersect(const Circle<T>& other) const
  {
    /*
     * See: Eric Hartmann - Geometry and Algorithms for COMPUTER AIDED DESIGN
     * https://www2.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
     */
    auto cp = c() - a() * other.center().x() - b() * other.center().y();
    auto discriminant = other.radius() * other.radius() * (
        a() * a() + b() * b()) - cp * cp;
    if (discriminant > T{0}) // Two intersections
    {
      return std::make_pair(
          other.center() +
          Point2D<T>(
            (a() * cp + b() * Rodin::Math::sqrt(discriminant)) / (a() * a() + b() *
              b()),
            (b() * cp - a() * Rodin::Math::sqrt(discriminant)) / (a() * a() + b() *
              b())
            ),
          other.center() +
          Point2D<T>(
            (a() * cp - b() * Rodin::Math::sqrt(discriminant)) / (a() * a() + b() *
              b()),
            (b() * cp + a() * Rodin::Math::sqrt(discriminant)) / (a() * a() + b() *
              b())
            )
          );
    }
    else if (discriminant == T{0}) // Tangent line
    {
      return other.center() + Point2D<T>(a() * cp, b() * cp);
    }
    else // No intersections
      return std::nullopt;
  }

  template <class T>
  inline
  constexpr
  std::optional<LineSegment2D<T>> Line2D<T>::connect(const Circle<T>& other) const
  {
    // Compute line segment from line to center
    auto ls = connect(other.center());

    if (ls)
    {
      if (other(ls->start()) > T{0}) // Closest point is outside the circle
      {
        // Line segment goes from the closest point on the line to the circle
        return LineSegment2D<T>(
            ls->start(),
            other.center() - other.radius() * ls->direction()
            );
      }
      else // Point is inside circle, line touches circle
      {
        return {};
      }
    }
    else // Line goes through center, touches the circle
    {
      return {};
    }
  }
}

#endif
