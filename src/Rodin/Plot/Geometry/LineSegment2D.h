/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_SEGMENT2D_H
#define RODIN_CORE_GEOMETRY_SEGMENT2D_H

#include <optional>

#include <Magnum/Math/Vector2.h>

#include "Rodin/Core/Common.h"

#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Plot::Geometry
{
  template <class T>
  class LineSegment2D : public Base<T, LineSegment2D<T>>
  {
   public:
    constexpr
    LineSegment2D(const Point2D<T>& start, const Point2D<T>& end);

    /**
     * Given \f$ t \in \mathbb{R} \f$, computes the parametric point
     * \f$ (x(t), y(t)) \f$ of the line segment parametrization:
     * \f[
     * (x(t), y(t)) = (x_0, y_0) + t \left((x_1, y_1) - (x_0, y_0)\right)
     * \f]
     * where
     * \f$(x_0, y_0)\f$ and \f$(x_1, y_1)\f$ are the start and end of the line
     * segment, respectively.
     *
     * @returns \f$ (x(t), x(t)) \f$
     */
    inline
    constexpr
    Point2D<T> operator()(const T& t) const;

    inline
    constexpr
    T length() const;

    /**
     * @returns The unit vector from `start()` to `end()`.
     */
    inline
    constexpr
    Magnum::Math::Vector2<T> direction() const;

    inline
    constexpr
    Point2D<T> start() const;

    inline
    constexpr
    Point2D<T> end() const;

    inline
    constexpr
    LineSegment2D& setStart(const Point2D<T>& start);

    inline
    constexpr
    LineSegment2D& setEnd(const Point2D<T>& end);

    inline
    constexpr
    LineSegment2D reverse() const;

    /**
     * Computes the intersection point (if any) of two line segments.
     * @retval std::nullopt If the line segments do not intersect or are
     * parallel with each other.
     * @retval Point2D If the line segments intersect
     */
    inline
    constexpr
    std::optional<Point2D<T>> intersect(const LineSegment2D<T>& other)
    const;

   private:
    Point2D<T>  m_start,
            m_end;
  };
}

#include "LineSegment2D.hpp"

#endif

