/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_EXAMPLES_WITNESS_BALL_H
#define RODIN_EXAMPLES_WITNESS_BALL_H

#include <Rodin/Math/SpatialVector.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_d.h>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace Examples
{
  /**
   * Checked smallest enclosing ball in reference dimensions one to three.
   *
   * Architecture: translate and normalize the input, call CGAL's smallest
   * enclosing sphere algorithm, then check containment and the convex-hull
   * optimality condition on active points. Degenerate failures use a finite
   * support search of at most dimension plus one points. All results are
   * floating-point solutions, not interval certificates.
   */
  class Ball
  {
    public:
      using Point = Rodin::Math::SpatialVector<Rodin::Real>;
      using Points = std::vector<Point>;
      using Real = Rodin::Real;
      static constexpr Real Tolerance = 1e-10;

      /** @param points Nonempty cell vertices. @param seed Input-order seed. */
      Ball(const Points& points, size_t seed)
        : m_dimension(points.front().size()), m_center(m_dimension)
      {
        Point origin(m_dimension);
        for (const auto& point : points)
          origin += point;
        origin /= Real(points.size());
        Real scale = 0;
        for (const auto& point : points)
          scale = std::max(scale, (point-origin).norm());
        if (scale == 0)
        {
          m_center = origin;
          return;
        }
        for (const auto& point : points)
          m_points.emplace_back((point-origin)/scale);
        switch (m_dimension)
        {
          case 1: enclose<1>(seed); break;
          case 2: enclose<2>(seed); break;
          case 3: enclose<3>(seed); break;
          default: throw std::runtime_error("Unsupported enclosing-ball dimension.");
        }
        if (!certifies())
        {
          m_fallback = true;
          if (!support())
            throw std::runtime_error("Smallest enclosing ball support is numerically unresolved.");
        }
        m_center = Point(origin+scale*m_center);
      }

      /** @return Checked enclosing-ball center. */
      const Point& getCenter() const { return m_center; }
      /** @return Whether a degenerate library result required support search. */
      bool usedFallback() const { return m_fallback; }

    private:
      template<int Dimension>
      void enclose(size_t seed)
      {
        using Kernel = CGAL::Cartesian_d<Real>;
        using Traits = CGAL::Min_sphere_of_spheres_d_traits_d<Kernel, Real, Dimension>;
        std::vector<typename Traits::Sphere> input;
        for (const auto& point : m_points)
        {
          const auto* coordinates = point.getData().data();
          input.emplace_back(typename Kernel::Point_d(Dimension, coordinates, coordinates+Dimension), 0);
        }
        std::mt19937_64 random(seed);
        std::shuffle(input.begin(), input.end(), random);
        CGAL::Min_sphere_of_spheres_d<Traits> ball(input.begin(), input.end());
        auto center = ball.center_cartesian_begin();
        for (size_t axis = 0; axis < m_dimension; ++axis)
          m_center(axis) = *center++;
        m_squaredRadius = ball.radius()*ball.radius();
      }

      bool contains() const
      {
        if (!std::isfinite(m_squaredRadius))
          return false;
        for (const auto& point : m_points)
        {
          if ((point-m_center).squaredNorm() > m_squaredRadius+Tolerance)
            return false;
        }
        return true;
      }

      bool certifies() const
      {
        if (!contains())
          return false;
        std::vector<size_t> active, selected;
        for (size_t i = 0; i < m_points.size(); ++i)
        {
          if (std::abs((m_points[i]-m_center).squaredNorm()-m_squaredRadius) <= Tolerance)
            active.push_back(i);
        }
        Eigen::VectorXd rhs(m_dimension+1);
        rhs.head(m_dimension) = m_center.getData().head(m_dimension);
        rhs(m_dimension) = 1;
        for (size_t count = 1; count <= m_dimension+1; ++count)
        {
          auto visit = [&](auto&& self, size_t begin) -> bool {
            if (selected.size() == count)
            {
              Eigen::MatrixXd matrix(m_dimension+1, count);
              for (size_t column = 0; column < count; ++column)
              {
                matrix.col(column).head(m_dimension) = m_points[selected[column]].getData().head(m_dimension);
                matrix(m_dimension, column) = 1;
              }
              const Eigen::VectorXd weights = matrix.completeOrthogonalDecomposition().solve(rhs);
              return weights.minCoeff() >= -Tolerance && (matrix*weights-rhs).norm() <= Tolerance;
            }
            for (size_t i = begin; i < active.size(); ++i)
            {
              selected.push_back(active[i]);
              if (self(self, i+1))
                return true;
              selected.pop_back();
            }
            return false;
          };
          if (visit(visit, 0))
            return true;
        }
        return false;
      }

      bool support()
      {
        std::vector<size_t> selected;
        for (size_t count = 1; count <= m_dimension+1; ++count)
        {
          auto visit = [&](auto&& self, size_t begin) -> bool {
            if (selected.size() == count)
            {
              const Point& origin = m_points[selected[0]];
              m_center = origin;
              if (count > 1)
              {
                Eigen::MatrixXd differences(count-1, m_dimension);
                Eigen::VectorXd rhs(count-1);
                for (size_t i = 1; i < count; ++i)
                {
                  differences.row(i-1) = (m_points[selected[i]]-origin).getData().head(m_dimension).transpose();
                  rhs(i-1) = differences.row(i-1).squaredNorm()/2;
                }
                const Eigen::MatrixXd gram = differences*differences.transpose();
                const Eigen::VectorXd weights = gram.completeOrthogonalDecomposition().solve(rhs);
                if (!weights.allFinite() || (gram*weights-rhs).norm() > Tolerance ||
                    weights.minCoeff() < -Tolerance || weights.sum() > 1+Tolerance)
                  return false;
                m_center = Point(origin.getData().head(m_dimension)+differences.transpose()*weights);
              }
              m_squaredRadius = (m_center-origin).squaredNorm();
              return contains();
            }
            for (size_t i = begin; i < m_points.size(); ++i)
            {
              selected.push_back(i);
              if (self(self, i+1))
                return true;
              selected.pop_back();
            }
            return false;
          };
          if (visit(visit, 0))
            return true;
        }
        return false;
      }

      size_t m_dimension;
      Point m_center;
      Real m_squaredRadius = 0;
      Points m_points;
      bool m_fallback = false;
  };
}
#endif
