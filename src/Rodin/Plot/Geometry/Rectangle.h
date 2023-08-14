/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_RECTANGLE_H
#define RODIN_CORE_GEOMETRY_RECTANGLE_H

#include "ForwardDecls.h"

namespace Rodin::Plot::Geometry
{
  template <class T>
  class Rectangle
  {
   public:
    enum Corner
    {
      BottomLeft,
      BottomRight,
      TopLeft,
      TopRight
    };

    constexpr
    Rectangle(const Point2D<T>& bottomLeft, const Point2D<T>& topRight);

    inline
    constexpr
    bool contains(const Point2D<T>& p) const;

    template <Corner c>
    inline
    constexpr
    Point2D<T> getCorner() const;

    inline
    constexpr
    T height() const;

    inline
    constexpr
    T width() const;

    inline
    constexpr
    T area() const;

   private:
    Point2D<T>  m_bottomLeft,
            m_topRight;
  };
}

#include "Rectangle.hpp"

#endif
