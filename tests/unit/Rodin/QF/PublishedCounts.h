/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Published point counts, and the counting bound they must respect.
 *
 * These are the targets the generators are measured against. They are recorded
 * once, here, because the two families publish *different* counts for the same
 * element and confusing them silently moves the goalposts: the Xiao--Gimbutas
 * tetrahedron of strength four has eleven points and no symmetry at all, while
 * the Witherden--Vincent tetrahedron of the same strength has fourteen and is
 * fully symmetric. Neither number is wrong, and a symmetric generator measured
 * against the asymmetric column looks like a failure when it is not.
 */
#ifndef RODIN_TESTS_UNIT_QF_PUBLISHEDCOUNTS_H
#define RODIN_TESTS_UNIT_QF_PUBLISHEDCOUNTS_H

#include <cstddef>
#include <map>
#include <vector>

#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Tests::QF
{
  /// @brief A published count, or absent where the family publishes none.
  using CountTable = std::map<size_t, size_t>;

  /**
   * @brief Builds a table from a degree-one-upwards list.
   * @param counts Counts for degrees 1, 2, 3, ... in order.
   */
  inline CountTable tabulate(const std::vector<size_t>& counts)
  {
    CountTable out;
    for (size_t i = 0; i < counts.size(); ++i)
      out[i + 1] = counts[i];
    return out;
  }

  /**
   * @brief Point counts of the fully symmetric, positive-weight, interior
   * rules of Witherden and Vincent.
   *
   * Table 1 of arXiv:1409.1865. The two-dimensional domains are published to
   * strength twenty, the three-dimensional ones to strength ten.
   */
  inline const CountTable& witherdenVincentCounts(Geometry::Polytope::Type g)
  {
    static const CountTable triangle = tabulate(
      {1, 3, 6, 6, 7, 12, 15, 16, 19, 25, 28, 33, 37, 42, 49, 55, 60, 67, 73, 79});
    static const CountTable quadrilateral = tabulate(
      {1, 4, 4, 8, 8, 12, 12, 20, 20, 28, 28, 37, 37, 48, 48, 60, 60, 72, 72, 85});
    static const CountTable tetrahedron = tabulate({1, 4, 8, 14, 14, 24, 35, 46, 59, 81});
    static const CountTable wedge = tabulate({1, 5, 8, 11, 16, 28, 35, 46, 60, 85});
    static const CountTable pyramid = tabulate({1, 5, 6, 10, 15, 24, 31, 47, 62, 83});
    static const CountTable hexahedron = tabulate({1, 6, 6, 14, 14, 34, 34, 58, 58, 90});
    static const CountTable none;

    switch (g)
    {
      case Geometry::Polytope::Type::Triangle:
        return triangle;
      case Geometry::Polytope::Type::Quadrilateral:
        return quadrilateral;
      case Geometry::Polytope::Type::Tetrahedron:
        return tetrahedron;
      case Geometry::Polytope::Type::Wedge:
        return wedge;
      case Geometry::Polytope::Type::Pyramid:
        return pyramid;
      case Geometry::Polytope::Type::Hexahedron:
        return hexahedron;
      default:
        return none;
    }
  }

  /**
   * @brief Point counts of the rules of Xiao and Gimbutas.
   *
   * Taken from the tables vendored by modepy, which reproduce
   * doi:10.1016/j.camwa.2009.10.027. These rules come from node elimination
   * and are asymmetric in general, so they are not a target a symmetry-based
   * generator can be expected to meet.
   *
   * The triangle is published to strength fifty and the tetrahedron only to
   * strength fifteen; beyond that the family offers no target at all, and a
   * rule is judged by the counting bound and by whichever rule is otherwise
   * best known.
   */
  inline const CountTable& xiaoGimbutasCounts(Geometry::Polytope::Type g)
  {
    static const CountTable triangle = tabulate(
      {1, 3, 6, 6, 7, 12, 15, 16, 19, 25, 28, 33, 37, 42, 49, 55, 60, 67, 73, 79});
    static const CountTable tetrahedron =
      tabulate({1, 4, 6, 11, 14, 23, 31, 44, 57, 74, 95, 122, 146, 177, 214});
    static const CountTable none;

    switch (g)
    {
      case Geometry::Polytope::Type::Triangle:
        return triangle;
      case Geometry::Polytope::Type::Tetrahedron:
        return tetrahedron;
      default:
        return none;
    }
  }

  /// @brief The published count, or zero where the family publishes none.
  inline size_t publishedCount(const CountTable& table, size_t degree)
  {
    const auto found = table.find(degree);
    return (found == table.end()) ? 0 : found->second;
  }

  /**
   * @brief Smallest point count any rule of this strength could have.
   *
   * A rule of strength @f$ p @f$ in @f$ d @f$ dimensions must reproduce
   * @f$ \binom{p+d}{d} @f$ independent moments, and each point contributes
   * @f$ d @f$ coordinates and one weight, so
   * @f[
   *   n \ge \left\lceil \binom{p+d}{d} \big/ (d+1) \right\rceil .
   * @f]
   * No rule can go below this, so a generated count that does is not a
   * discovery but a bug --- most often a moment system that has stopped
   * spanning the space it is supposed to.
   */
  inline size_t countingBound(size_t dimension, size_t degree)
  {
    size_t moments = 1;
    for (size_t i = 1; i <= dimension; ++i)
      moments = moments * (degree + i) / i;
    return (moments + dimension) / (dimension + 1);
  }
}

#endif
