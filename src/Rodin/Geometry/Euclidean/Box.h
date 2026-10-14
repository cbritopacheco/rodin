/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_BOX_H
#define RODIN_GEOMETRY_EUCLIDEAN_BOX_H

#include <array>
#include <algorithm>
#include "Point.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Closed axis-aligned box @f$ \prod_i[l_i,u_i] @f$.
   * @tparam Dimension Nonnegative ambient dimension; the zero-dimensional box has counting measure one.
   * @tparam T Real coordinate scalar type.
   */
  template <size_t Dimension, class T>
  class Box : public Base<Box<Dimension, T>, T>
  {
    public:
      /// @brief Planar corner names retained for the rectangle alias.
      enum Corner
      {
        BottomLeft,
        BottomRight,
        TopLeft,
        TopRight
      };
      /// @brief Constructs a box. @param lower Componentwise lower corner. @param upper Componentwise upper corner.
      Box(const Point<Dimension, T>& lower, const Point<Dimension, T>& upper)
        : m_lower(lower),
          m_upper(upper)
      {
        for (size_t i = 0; i < Dimension; ++i)
        {
          if (!(lower[i] <= upper[i]))
            Alert::MemberFunctionException(*this, __func__)
              << "Box requires ordered bounds in coordinate " << i << Alert::Raise;
        }
      }
      /// @brief Gets the ambient dimension. @returns Coordinate count.
      static constexpr size_t getDimension()
      {
        return Dimension;
      }
      /// @brief Gets the lower corner. @returns Componentwise minima.
      const Point<Dimension, T>& lower() const
      {
        return m_lower;
      }
      /// @brief Gets the upper corner. @returns Componentwise maxima.
      const Point<Dimension, T>& upper() const
      {
        return m_upper;
      }
      /// @brief Tests closed containment. @param point Point. @returns Whether all coordinates are in bounds.
      bool contains(const Point<Dimension, T>& point) const
      {
        for (size_t i = 0; i < Dimension; ++i)
        {
          if (!(m_lower[i] <= point[i] && point[i] <= m_upper[i]))
            return false;
        }
        return true;
      }
      /// @brief Selects a corner. @param upperCoordinates True selects the upper bound in that coordinate. @returns Corner point.
      Point<Dimension, T> getCorner(
        const std::array<bool, Dimension>& upperCoordinates) const
      {
        auto point = m_lower;
        for (size_t i = 0; i < Dimension; ++i)
        {
          if (upperCoordinates[i])
            point[i] = m_upper[i];
        }
        return point;
      }
      /// @brief Selects a planar corner. @tparam C Corner identifier. @returns Corner point.
      template <Corner C>
      Point<Dimension, T> getCorner() const
        requires(Dimension == 2)
      {
        static_assert(C >= BottomLeft && C <= TopRight);
        return getCorner(
          {C == BottomRight || C == TopRight, C == TopLeft || C == TopRight});
      }
      /// @brief Computes a side length. @param coordinate Coordinate index. @returns Extent in that direction.
      T extent(size_t coordinate) const
      {
        if (coordinate >= Dimension)
          Alert::MemberFunctionException(*this, __func__)
            << "Box coordinate is outside its dimension" << Alert::Raise;
        return m_upper[coordinate] - m_lower[coordinate];
      }
      /// @brief Computes measure. @returns Product of all side lengths.
      T measure() const
      {
        if constexpr (Dimension == 0)
          return T(1);
        else if constexpr (Dimension == 1)
          return extent(0);
        else if constexpr (Dimension == 2)
          return extent(0) * extent(1);
        else if constexpr (Dimension == 3)
          return extent(0) * extent(1) * extent(2);
        else
        {
          T result = T(1);
          for (size_t i = 0; i < Dimension; ++i)
            result *= extent(i);
          return result;
        }
      }
      /// @brief Projects onto the closed box. @param point Point. @returns Closest point in the box.
      Point<Dimension, T> project(const Point<Dimension, T>& point) const
      {
        auto result = point;
        for (size_t i = 0; i < Dimension; ++i)
          result[i] = std::clamp(point[i], m_lower[i], m_upper[i]);
        return result;
      }
      /// @brief Computes distance to the box. @param point Point. @returns Euclidean distance, zero inside.
      T distance(const Point<Dimension, T>& point) const
      {
        return (point - project(point)).norm();
      }
      /// @brief Computes box intersection. @param other Box. @returns Closed intersection, including degenerate touching boxes.
      Optional<Box> intersect(const Box& other) const
      {
        auto lower = m_lower, upper = m_upper;
        for (size_t i = 0; i < Dimension; ++i)
        {
          lower[i] = std::max(lower[i], other.m_lower[i]);
          upper[i] = std::min(upper[i], other.m_upper[i]);
          if (lower[i] > upper[i])
            return {};
        }
        return Box(lower, upper);
      }
      /// @brief Gets planar width. @returns First side length.
      T width() const
        requires(Dimension == 2)
      {
        return extent(0);
      }
      /// @brief Gets planar height. @returns Second side length.
      T height() const
        requires(Dimension == 2)
      {
        return extent(1);
      }
      /// @brief Gets planar area. @returns Product of both side lengths.
      T area() const
        requires(Dimension == 2)
      {
        return measure();
      }

    private:
      Point<Dimension, T> m_lower, m_upper;
  };
}

#endif
