/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_POINT2D_H
#define RODIN_GEOMETRY_EUCLIDEAN_POINT2D_H

/**
 * @file
 * @brief Point in 2D Euclidean space.
 */

#include <initializer_list>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Math/SpatialVector.h"
#include "Base.h"

#include "Rodin/Types.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Point in 2D Euclidean space.
   *
   * Represents a point @f$ p = (x, y) @f$ in two-dimensional Euclidean space.
   * Coordinates use the fixed-capacity Math::SpatialVector storage with two
   * active components.
   *
   * @tparam T Scalar type (e.g., float, double)
   *
   * # Operations
   *
   * The class supports:
   * - Spatial displacements, coordinate addition and scaling
   * - Connection to other geometric objects (Line2D, Circle, other points)
   * - Distance computations
   *
   *
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line2_d.html">Line2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_circle.html">Circle</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
   */
  template <class T>
  class Point2D : public Base<Point2D<T>, T>
  {
    public:
      /// @brief Spatial coordinate storage type.
      using Vector = Math::SpatialVector<T>;

      /// @brief Constructs the origin with two active coordinates.
      constexpr Point2D()
        : m_coordinates(2)
      {}

      /// @brief Constructs a point. @param x First coordinate. @param y Second coordinate.
      constexpr Point2D(T x, T y)
        : m_coordinates({x, y})
      {}

      /// @brief Constructs a point from two coordinates. @param coordinates Initial coordinates.
      constexpr Point2D(std::initializer_list<T> coordinates)
        : Point2D()
      {
        requireTwoCoordinates(coordinates.size());
        x() = *coordinates.begin();
        y() = *(coordinates.begin() + 1);
      }

      /// @brief Constructs a point from a spatial vector. @param coordinates Two active coordinates.
      constexpr Point2D(const Vector& coordinates)
        : Point2D()
      {
        *this = coordinates;
      }

      /// @brief Assigns two spatial coordinates. @param coordinates New coordinates. @returns This point.
      Point2D& operator=(const Vector& coordinates)
      {
        requireTwoCoordinates(coordinates.size());
        m_coordinates = coordinates;
        return *this;
      }

      /// @brief Gets the mutable spatial coordinates. @returns Storage with two active components; callers must preserve its size.
      Vector& asVector()
      {
        return m_coordinates;
      }

      /// @brief Gets the spatial coordinates. @returns Storage with two active components.
      const Vector& asVector() const
      {
        return m_coordinates;
      }

      /// @brief Gets the mutable first coordinate. @returns The first coordinate.
      T& x()
      {
        return m_coordinates.x();
      }
      /// @brief Gets the first coordinate. @returns The first coordinate.
      const T& x() const
      {
        return m_coordinates.x();
      }
      /// @brief Gets the mutable second coordinate. @returns The second coordinate.
      T& y()
      {
        return m_coordinates.y();
      }
      /// @brief Gets the second coordinate. @returns The second coordinate.
      const T& y() const
      {
        return m_coordinates.y();
      }

      /// @brief Gets a mutable coordinate. @param i Coordinate index. @returns The coordinate.
      T& operator[](size_t i)
      {
        assert(i < 2);
        return m_coordinates[i];
      }
      /// @brief Gets a coordinate. @param i Coordinate index. @returns The coordinate.
      const T& operator[](size_t i) const
      {
        assert(i < 2);
        return m_coordinates[i];
      }

      /// @brief Gets the distance from the origin. @returns Euclidean coordinate norm.
      auto norm() const
      {
        return m_coordinates.norm();
      }

      /// @brief Compares coordinates. @param other Other point. @returns Whether both coordinates are equal.
      bool operator==(const Point2D& other) const
      {
        return x() == other.x() && y() == other.y();
      }

      /// @brief Adds coordinates. @param lhs First point. @param rhs Second point. @returns Coordinate sum.
      friend Point2D operator+(const Point2D& lhs, const Point2D& rhs)
      {
        return Point2D(lhs.x() + rhs.x(), lhs.y() + rhs.y());
      }

      /// @brief Computes a displacement. @param lhs End point. @param rhs Start point. @returns Spatial displacement.
      friend Vector operator-(const Point2D& lhs, const Point2D& rhs)
      {
        return Vector{lhs.x() - rhs.x(), lhs.y() - rhs.y()};
      }

      /// @brief Translates a point. @param point Initial point. @param offset Spatial displacement. @returns Translated point.
      friend Point2D operator+(const Point2D& point, const Vector& offset)
      {
        point.requireTwoCoordinates(offset.size());
        return Point2D(point.x() + offset.x(), point.y() + offset.y());
      }

      /// @brief Translates a point. @param offset Spatial displacement. @param point Initial point. @returns Translated point.
      friend Point2D operator+(const Vector& offset, const Point2D& point)
      {
        return point + offset;
      }

      /// @brief Subtracts a displacement. @param point Initial point. @param offset Spatial displacement. @returns Translated point.
      friend Point2D operator-(const Point2D& point, const Vector& offset)
      {
        point.requireTwoCoordinates(offset.size());
        return Point2D(point.x() - offset.x(), point.y() - offset.y());
      }

      /// @brief Scales coordinates. @param point Point to scale. @param scalar Scale factor. @returns Scaled point.
      friend Point2D operator*(const Point2D& point, const T& scalar)
      {
        return Point2D(point.x() * scalar, point.y() * scalar);
      }

      /// @brief Scales coordinates. @param scalar Scale factor. @param point Point to scale. @returns Scaled point.
      friend Point2D operator*(const T& scalar, const Point2D& point)
      {
        return point * scalar;
      }

      /// @brief Divides coordinates. @param point Point to scale. @param scalar Nonzero divisor. @returns Scaled point.
      friend Point2D operator/(const Point2D& point, const T& scalar)
      {
        return Point2D(point.x() / scalar, point.y() / scalar);
      }

      /**
       * @brief Computes the line segment connecting this point to another.
       *
       * Creates a line segment where the start point is the current point
       * and the end point is the other point.
       *
       * @param[in] other Point to connect to
       * @returns Optional line segment connecting the points
       * @retval std::nullopt If the points are equal (no segment can be formed)
       * @retval LineSegment2D If the points are distinct
       *
       * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
       */
      inline
      constexpr
      Optional<LineSegment2D<T>> connect(const Point2D& other) const;

      /**
       * @brief Computes the line segment connecting this point to nearest point on a line.
       *
       * Creates a line segment where the start point is the current point
       * and the end point is the nearest point on the given line.
       *
       * @param[in] other Line to connect to
       * @returns Optional line segment
       * @retval std::nullopt If the point lies on the line
       * @retval LineSegment2D Perpendicular segment from point to line
       *
       * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line2_d.html">Line2D</a>
       * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
       */
      inline
      constexpr
      Optional<LineSegment2D<T>> connect(const Line2D<T>& other) const;

      /**
       * @brief Computes the line segment connecting this point to nearest point on a circle.
       *
       * Creates a line segment where the start point is the current point and
       * the end point is the nearest point on the circle's circumference.
       *
       * @param[in] other Circle to connect to
       * @returns Optional line segment
       * @retval std::nullopt If the point lies on the circle
       * @retval LineSegment2D Segment from point to nearest point on circle
       *
       * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_circle.html">Circle</a>
       * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
       */
      inline
      constexpr
      Optional<LineSegment2D<T>> connect(const Circle<T>& other) const;

    private:
      /// @brief Validates the coordinate dimension. @param size Number of active coordinates.
      void requireTwoCoordinates(size_t size) const
      {
        if (size != 2)
          Alert::MemberFunctionException(*this, __func__)
            << "Point2D requires two coordinates; received " << size << Alert::Raise;
      }

      Vector m_coordinates;
  };
}

#include "Point2D.hpp"

#endif
