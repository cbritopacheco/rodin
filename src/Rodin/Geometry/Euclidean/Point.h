/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_POINT_H
#define RODIN_GEOMETRY_EUCLIDEAN_POINT_H

/**
 * @file
 * @brief Point with a compile-time Euclidean dimension.
 */

#include <initializer_list>
#include <limits>
#include <cmath>
#include <array>
#include <type_traits>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Math/Vector.h"
#include "Base.h"

#include "Rodin/Types.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Point @f$ p\in\mathbb{R}^{Dimension} @f$ with compile-time fixed coordinate storage.
   * Arithmetic applies componentwise; subtracting two points yields a spatial
   * displacement. The template dimension is a nonnegative compile-time coordinate count.
   * @tparam Dimension Ambient coordinate count.
   * @tparam T Coordinate scalar type.
   */
  template <size_t Dimension, class T>
  class Point : public Base<Point<Dimension, T>, T>
  {
    public:
      /// @brief Compile-time coordinate storage type.
      using Vector = Math::FixedSizeVector<T, Dimension>;
      static_assert(Dimension <= static_cast<size_t>(std::numeric_limits<int>::max()));

      /// @brief Gets the compile-time coordinate dimension. @returns The active coordinate count.
      static constexpr size_t getDimension() noexcept
      {
        return Dimension;
      }

      /// @brief Constructs the origin with Dimension active coordinates.
      constexpr Point()
        : m_coordinates(Vector::Zero())
      {}

      /// @brief Constructs a point from its coordinates. @param coordinates Coordinate values.
      template <class... Coordinates>
        requires(sizeof...(Coordinates) == Dimension &&
          (std::is_convertible_v<Coordinates, T> && ...))
      constexpr Point(Coordinates... coordinates)
        : Point()
      {
        const std::array<T, Dimension> values{static_cast<T>(coordinates)...};
        for (size_t i = 0; i < Dimension; ++i)
          m_coordinates[i] = values[i];
      }
      /// @brief Constructs a point from Dimension coordinates. @param coordinates Initial coordinates.
      constexpr Point(std::initializer_list<T> coordinates)
        : Point()
      {
        require(coordinates.size());
        size_t i = 0;
        for (const T& coordinate : coordinates)
          m_coordinates[i++] = coordinate;
      }
      /// @brief Copies a column-vector expression into fixed storage. @param coordinates Dimension coordinate values.
      template <class Derived>
        requires((Derived::RowsAtCompileTime == Eigen::Dynamic ||
                   Derived::RowsAtCompileTime == Dimension) &&
          (Derived::ColsAtCompileTime == Eigen::Dynamic ||
            Derived::ColsAtCompileTime == 1))
      explicit Point(const Eigen::MatrixBase<Derived>& coordinates)
        : Point()
      {
        *this = coordinates;
      }
      /// @brief Assigns a column-vector expression. @param coordinates New coordinates. @returns This point.
      template <class Derived>
        requires((Derived::RowsAtCompileTime == Eigen::Dynamic ||
                   Derived::RowsAtCompileTime == Dimension) &&
          (Derived::ColsAtCompileTime == Eigen::Dynamic ||
            Derived::ColsAtCompileTime == 1))
      Point& operator=(const Eigen::MatrixBase<Derived>& coordinates)
      {
        require(coordinates.size());
        if (coordinates.cols() != 1)
          Alert::MemberFunctionException(*this, __func__)
            << "Point requires a column vector" << Alert::Raise;
        const Vector owned = coordinates.template cast<T>();
        m_coordinates = owned;
        return *this;
      }

      /// @brief Gets the mutable fixed-size coordinates. @returns Fixed storage with Dimension components.
      Vector& asVector()
      {
        return m_coordinates;
      }

      /// @brief Gets the fixed-size coordinates. @returns Fixed storage with Dimension components.
      const Vector& asVector() const
      {
        return m_coordinates;
      }

      /// @brief Gets the mutable first coordinate. @returns The first coordinate.
      T& x()
        requires(Dimension >= 1)
      {
        return m_coordinates.x();
      }
      /// @brief Gets the first coordinate. @returns The first coordinate.
      const T& x() const
        requires(Dimension >= 1)
      {
        return m_coordinates.x();
      }
      /// @brief Gets the mutable second coordinate. @returns The second coordinate.
      T& y()
        requires(Dimension >= 2)
      {
        return m_coordinates.y();
      }
      /// @brief Gets the second coordinate. @returns The second coordinate.
      const T& y() const
        requires(Dimension >= 2)
      {
        return m_coordinates.y();
      }

      /// @brief Gets the mutable third coordinate. @returns The third coordinate.
      T& z()
        requires(Dimension >= 3)
      {
        return m_coordinates.z();
      }
      /// @brief Gets the third coordinate. @returns The third coordinate.
      const T& z() const
        requires(Dimension >= 3)
      {
        return m_coordinates.z();
      }

      /// @brief Gets a mutable coordinate. @param i Coordinate index. @returns The coordinate.
      T& operator[](size_t i)
      {
        assert(i < Dimension);
        return m_coordinates[i];
      }
      /// @brief Gets a coordinate. @param i Coordinate index. @returns The coordinate.
      const T& operator[](size_t i) const
      {
        assert(i < Dimension);
        return m_coordinates[i];
      }

      /// @brief Gets the distance from the origin. @returns Euclidean coordinate norm.
      auto norm() const
      {
        if constexpr (Dimension == 0)
          return T(0);
        else if constexpr (Dimension == 1)
          return std::abs(m_coordinates[0]);
        else
          return m_coordinates.norm();
      }

      /// @brief Compares coordinates. @param other Other point. @returns Whether both coordinates are equal.
      bool operator==(const Point& other) const
      {
        for (size_t i = 0; i < Dimension; ++i)
        {
          if (m_coordinates[i] != other.m_coordinates[i])
            return false;
        }
        return true;
      }

      /// @brief Adds coordinates. @param lhs First point. @param rhs Second point. @returns Coordinate sum.
      friend Point operator+(const Point& lhs, const Point& rhs)
      {
        return Point(lhs.m_coordinates + rhs.m_coordinates);
      }

      /// @brief Computes a displacement. @param lhs End point. @param rhs Start point. @returns Spatial displacement.
      friend Vector operator-(const Point& lhs, const Point& rhs)
      {
        return lhs.m_coordinates - rhs.m_coordinates;
      }

      /// @brief Translates a point. @param point Initial point. @param offset Column-vector displacement. @returns Translated point.
      template <class Derived>
      friend Point operator+(const Point& point, const Eigen::MatrixBase<Derived>& offset)
      {
        const Point displacement(offset);
        return Point(point.m_coordinates + displacement.m_coordinates);
      }
      /// @brief Translates a point. @param offset Column-vector displacement. @param point Initial point. @returns Translated point.
      template <class Derived>
      friend Point operator+(const Eigen::MatrixBase<Derived>& offset, const Point& point)
      {
        return point + offset;
      }
      /// @brief Subtracts a displacement. @param point Initial point. @param offset Column-vector displacement. @returns Translated point.
      template <class Derived>
      friend Point operator-(const Point& point, const Eigen::MatrixBase<Derived>& offset)
      {
        const Point displacement(offset);
        return Point(point.m_coordinates - displacement.m_coordinates);
      }

      /// @brief Scales coordinates. @param point Point to scale. @param scalar Scale factor. @returns Scaled point.
      friend Point operator*(const Point& point, const T& scalar)
      {
        return Point(point.m_coordinates * scalar);
      }

      /// @brief Scales coordinates. @param scalar Scale factor. @param point Point to scale. @returns Scaled point.
      friend Point operator*(const T& scalar, const Point& point)
      {
        return point * scalar;
      }

      /// @brief Divides coordinates. @param point Point to scale. @param scalar Nonzero divisor. @returns Scaled point.
      friend Point operator/(const Point& point, const T& scalar)
      {
        return Point(point.m_coordinates / scalar);
      }

      /// @brief Connects two distinct points. @param other End point. @returns Segment, absent for equal points.
      inline constexpr Optional<LineSegment<Dimension, T>> connect(
        const Point& other) const;
      /// @brief Connects to the closest line point. @param other Line. @returns Perpendicular segment, absent for a point on the line.
      inline constexpr Optional<LineSegment<Dimension, T>> connect(
        const Line<Dimension, T>& other) const;
      /// @brief Connects to the nearest sphere point. @param other Sphere. @returns Segment, absent at the center or on the surface where no unique nonzero connection exists.
      inline constexpr Optional<LineSegment<Dimension, T>> connect(
        const Sphere<Dimension, T>& other) const;

    private:
      /// @brief Validates the coordinate dimension. @param size Number of active coordinates.
      void require(size_t size) const
      {
        if (size != Dimension)
          Alert::MemberFunctionException(*this, __func__)
            << "Point requires " << Dimension << " coordinates; received " << size
            << Alert::Raise;
      }

      Vector m_coordinates;
  };
}

#include "Point.hpp"

#endif
