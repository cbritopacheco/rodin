/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_COMMON_H
#define RODIN_GEOMETRY_EUCLIDEAN_COMMON_H

/**
 * @file
 * @brief Common utility functions for Euclidean geometry.
 */

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Computes the barycenter (centroid) of multiple points.
   *
   * Given points @f$ p_0, p_1, \ldots, p_n @f$, computes:
   * @f[
   *   \text{barycenter} = \frac{1}{n+1} \sum_{i=0}^{n} p_i
   * @f]
   *
   * @tparam P0 Type of first point
   * @tparam Ps Types of remaining points
   * @param[in] p0 First point
   * @param[in] ps Remaining points
   * @returns Barycenter of all points
   *
   * @note All point types must support addition and scalar division.
   *
   * @par Example
   * @code{.cpp}
   * Point2D<double> p1{0, 0}, p2{2, 0}, p3{1, 2};
   * auto center = barycenter(p1, p2, p3); // Returns {1, 2/3}
   * @endcode
   */
  template <class P0, class ... Ps>
  auto barycenter(const P0& p0, const Ps&... ps)
  {
    return (p0 + ... + ps) / (sizeof...(ps) + 1);
  }
}

#endif
