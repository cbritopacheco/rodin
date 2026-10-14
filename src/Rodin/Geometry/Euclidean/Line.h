/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_LINE_H
#define RODIN_GEOMETRY_EUCLIDEAN_LINE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <variant>
#include "Rodin/Math/Matrix.h"
#include "Point.h"
#include "LineSegment.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Affine line in a compile-time Euclidean dimension.
   * @f$ L(t)=p+tv @f$ with nonzero direction @f$ v @f$.
   * Projection uses the Euclidean inner product. Closest points use direct formulas
   * through dimension three and a fixed-size QR solve in higher dimensions.
   * Unit direction columns use a rank threshold of twice scalar machine epsilon
   * (relative to the largest QR pivot in higher dimensions). Numerically parallel directions
   * have rank below two and admit no unique line intersection.
   * Intersections return only unique points; coincident lines return no point.
   * @tparam Dimension Positive ambient dimension.
   * @tparam T Real coordinate scalar type.
   */
  template <size_t Dimension, class T>
  class Line : public Base<Line<Dimension, T>, T>
  {
    public:
      using Vector = typename Point<Dimension, T>::Vector;
      using Intersections = std::variant<std::nullopt_t, Point<Dimension, T>,
        std::pair<Point<Dimension, T>, Point<Dimension, T>>>;

      /// @brief Constructs a line through two points. @param p First point. @param q Distinct second point.
      Line(const Point<Dimension, T>& p, const Point<Dimension, T>& q)
        : Line(p, q - p)
      {}

      /// @brief Constructs an affine line. @param origin Origin. @param direction Nonzero direction.
      template <class Derived>
      Line(const Point<Dimension, T>& origin, const Eigen::MatrixBase<Derived>& direction)
        : m_origin(origin),
          m_direction(Vector::Zero()),
          m_coefficients{}
      {
        if (direction.size() != Dimension || direction.cols() != 1)
          Alert::MemberFunctionException(*this, __func__)
            << "Line requires a direction with " << Dimension << " coordinates"
            << Alert::Raise;
        m_direction = direction.template cast<T>();
        if (!(m_direction.array() != T(0)).any())
          Alert::MemberFunctionException(*this, __func__)
            << "Line requires a nonzero direction" << Alert::Raise;
        if constexpr (Dimension == 2)
          m_coefficients = {-m_direction[1], m_direction[0],
            -m_direction[1] * origin[0] + m_direction[0] * origin[1]};
      }

      /// @brief Constructs a supporting line. @param segment Nondegenerate segment.
      Line(const LineSegment<Dimension, T>& segment)
        : Line(segment.start(), segment.end())
      {}

      /// @brief Constructs a planar line @f$ ax+by=c @f$. @param a First coefficient. @param b Second coefficient. @param c Constant.
      Line(T a, T b, T c)
        requires(Dimension == 2)
        : m_origin(),
          m_direction(b, -a),
          m_coefficients{a, b, c}
      {
        if (a == T(0) && b == T(0))
          Alert::MemberFunctionException(*this, __func__)
            << "Line requires a nonzero normal" << Alert::Raise;
        if (std::abs(a) >= std::abs(b))
          m_origin[0] = c / a;
        else
          m_origin[1] = c / b;
      }

      /// @brief Constructs a planar line from three coefficients. @param coefficients Ordered coefficients a, b, c.
      Line(std::initializer_list<T> coefficients)
        requires(Dimension == 2)
        : m_origin(),
          m_direction(Vector::Zero()),
          m_coefficients{}
      {
        if (coefficients.size() != 3)
          Alert::MemberFunctionException(*this, __func__)
            << "Line requires three coefficients" << Alert::Raise;
        auto it = coefficients.begin();
        const T a = *it++, b = *it++, c = *it;
        *this = Line(a, b, c);
      }

      /// @brief Gets the ambient dimension. @returns Coordinate count.
      static constexpr size_t getDimension()
      {
        return Dimension;
      }
      /// @brief Gets the origin. @returns An anchor point.
      const Point<Dimension, T>& origin() const
      {
        return m_origin;
      }
      /// @brief Gets a unit direction. @returns Normalized direction.
      Vector direction() const
      {
        if constexpr (Dimension == 0)
          return Vector::Zero();
        else if constexpr (Dimension == 1)
        {
          Vector result;
          result[0] = m_direction[0] > T(0) ? T(1) : T(-1);
          return result;
        }
        else
          return m_direction / m_direction.stableNorm();
      }
      /// @brief Evaluates the affine parameterization. @param t Parameter. @returns Point on the line.
      Point<Dimension, T> pointAt(T t) const
      {
        return m_origin + t * m_direction;
      }
      /// @brief Projects a point orthogonally. @param p Point. @returns Closest point on the line.
      Point<Dimension, T> project(const Point<Dimension, T>& p) const
      {
        if constexpr (Dimension <= 1)
          return p;
        else
        {
          const Vector unit = direction();
          return m_origin + unit.dot(p - m_origin) * unit;
        }
      }
      /// @brief Computes point distance. @param p Point. @returns Euclidean distance.
      T distance(const Point<Dimension, T>& p) const
      {
        return (p - project(p)).norm();
      }

      /// @brief Computes the closest point on this line to another line. @param other Line. @returns Closest point; the origin is chosen for numerically parallel lines.
      Point<Dimension, T> closestPoint(const Line& other) const
      {
        const auto point = getUniqueClosestPoint(other);
        return point ? *point : m_origin;
      }
      /// @brief Computes line distance, including skew lines. @param other Line. @returns Minimum distance.
      T distance(const Line& other) const
      {
        if constexpr (Dimension <= 1)
          return T(0);
        else if constexpr (Dimension == 2)
        {
          if (getUniqueClosestPoint(other))
            return T(0);
          return other.distance(m_origin);
        }
        else if constexpr (Dimension == 3)
        {
          const Vector normal = direction().cross(other.direction());
          const T length = normal.norm();
          if (length <= rankTolerance())
            return other.distance(m_origin);
          return std::abs((other.m_origin - m_origin).dot(normal)) / length;
        }
        else
          return other.distance(closestPoint(other));
      }

      /**
       * @brief Computes a unique line intersection.
       * A candidate is accepted when the projection residual is at most
       * @f$ Dimension\,\epsilon\max(1,\|p\|,\|q\|,\|x\|) @f$.
       * This first-order rounding tolerance is a numerical policy, not an exact predicate.
       * @param other Other line. @returns Unique intersection, or no point for parallel or skew lines.
       */
      Optional<Point<Dimension, T>> intersect(const Line& other) const
      {
        const auto candidate = getUniqueClosestPoint(other);
        if (!candidate)
          return {};
        const auto& point = *candidate;
        const T scale =
          std::max({T(1), m_origin.norm(), other.m_origin.norm(), point.norm()});
        const T tolerance = T(Dimension) * std::numeric_limits<T>::epsilon() * scale;
        if (other.distance(point) > tolerance)
          return {};
        return point;
      }
      /// @brief Intersects a sphere. @param sphere Sphere. @returns Zero, one, or two points.
      Intersections intersect(const Sphere<Dimension, T>& sphere) const;
      /// @brief Connects to a disjoint sphere. @param sphere Sphere. @returns Unique shortest segment when separated.
      Optional<LineSegment<Dimension, T>> connect(
        const Sphere<Dimension, T>& sphere) const;

      /// @brief Gets the first planar normal coefficient. @returns Coefficient a.
      T a() const
        requires(Dimension == 2)
      {
        return m_coefficients[0];
      }
      /// @brief Gets the second planar normal coefficient. @returns Coefficient b.
      T b() const
        requires(Dimension == 2)
      {
        return m_coefficients[1];
      }
      /// @brief Gets the planar constant. @returns Coefficient c.
      T c() const
        requires(Dimension == 2)
      {
        return m_coefficients[2];
      }
      /// @brief Evaluates the planar implicit equation. @param p Point. @returns Signed equation value.
      T operator()(const Point<Dimension, T>& p) const
        requires(Dimension == 2)
      {
        return a() * p[0] + b() * p[1] - c();
      }
      /// @brief Evaluates planar ordinate. @param x Abscissa. @returns Ordinate, absent for vertical lines.
      Optional<T> operator()(T x) const
        requires(Dimension == 2)
      {
        if (b() == T(0))
          return {};
        return (c() - a() * x) / b();
      }
      /// @brief Gets planar slope. @returns Slope, absent for vertical lines.
      Optional<T> slope() const
        requires(Dimension == 2)
      {
        if (b() == T(0))
          return {};
        return -a() / b();
      }
      /// @brief Gets planar abscissa intercept. @returns Intercept when unique.
      Optional<T> xIntercept() const
        requires(Dimension == 2)
      {
        if (a() == T(0))
          return {};
        return c() / a();
      }
      /// @brief Gets planar ordinate intercept. @returns Intercept when unique.
      Optional<T> yIntercept() const
        requires(Dimension == 2)
      {
        if (b() == T(0))
          return {};
        return c() / b();
      }
      /// @brief Tests planar verticality. @returns Whether the second coefficient vanishes.
      bool isVertical() const
        requires(Dimension == 2)
      {
        return b() == T(0);
      }
      /// @brief Tests planar horizontality. @returns Whether the first coefficient vanishes.
      bool isHorizontal() const
        requires(Dimension == 2)
      {
        return a() == T(0);
      }

    private:
      /// @brief Gets the normalized-direction rank tolerance. @returns Twice scalar machine epsilon, matching the two-column QR pivot threshold.
      static constexpr T rankTolerance()
      {
        return T(2) * std::numeric_limits<T>::epsilon();
      }

      /**
       * @brief Solves the unique closest-point problem for two nonparallel lines.
       * Unit direction columns make the rank decision invariant under line
       * parameter rescaling. Dimensions zero and one have no rank-two span;
       * dimension two uses a determinant, dimension three a cross product,
       * and higher dimensions column-pivoted Householder QR on fixed storage.
       * @param other Other line.
       * @returns Unique closest point on this line, absent for numerically parallel directions.
       */
      Optional<Point<Dimension, T>> getUniqueClosestPoint(const Line& other) const
      {
        if constexpr (Dimension <= 1)
          return {};
        else
        {
          const Vector u = direction(), v = other.direction();
          const Vector offset = other.m_origin - m_origin;
          if constexpr (Dimension == 2)
          {
            const T determinant = u[0] * v[1] - u[1] * v[0];
            if (std::abs(determinant) <= rankTolerance())
              return {};
            const T parameter = (offset[0] * v[1] - offset[1] * v[0]) / determinant;
            return m_origin + parameter * u;
          }
          else if constexpr (Dimension == 3)
          {
            const Vector normal = u.cross(v);
            if (normal.norm() <= rankTolerance())
              return {};
            const T parameter = offset.cross(v).dot(normal) / normal.squaredNorm();
            return m_origin + parameter * u;
          }
          else
          {
            Math::FixedSizeMatrix<T, Dimension, 2> directions;
            directions.col(0) = u;
            directions.col(1) = -v;
            const Eigen::ColPivHouseholderQR<decltype(directions)> factorization(
              directions);
            if (factorization.rank() < 2)
              return {};
            const Math::FixedSizeVector<T, 2> parameters = factorization.solve(offset);
            return m_origin + parameters[0] * u;
          }
        }
      }

      Point<Dimension, T> m_origin;
      Vector m_direction;
      std::array<T, Dimension == 2 ? 3 : 0> m_coefficients;
  };
}

#include "Line.hpp"

#endif
