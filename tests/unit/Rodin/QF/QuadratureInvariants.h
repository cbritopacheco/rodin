/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_UNIT_QF_QUADRATUREINVARIANTS_H
#define RODIN_TESTS_UNIT_QF_QUADRATUREINVARIANTS_H

/**
 * @file
 * @brief Reference moments and rule invariants shared by the quadrature tests.
 *
 * The exactness check compares a rule against closed-form moments of the
 * reference polytope, so it needs no reference implementation and cannot
 * inherit an error from one. Every reference element below is the unit
 * element used by Rodin::QF; the moment formulae are derived for exactly
 * those domains.
 */

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Tests::QF
{
  /// @brief @f$ n! @f$ as a Real, for the small @f$ n @f$ used by moments.
  inline Real factorial(size_t n)
  {
    Real r = 1;
    for (size_t i = 2; i <= n; ++i)
      r *= static_cast<Real>(i);
    return r;
  }

  /**
   * @brief Exact value of @f$ \int_K x^a y^b z^c \, dK @f$ on the reference
   * polytope @p g.
   *
   * Domains, matching Rodin::QF's unit reference elements:
   * - Segment @f$ [0,1] @f$, Quadrilateral @f$ [0,1]^2 @f$,
   *   Hexahedron @f$ [0,1]^3 @f$: separable, @f$ \prod 1/(e_i+1) @f$.
   * - Triangle, Tetrahedron (unit simplices):
   *   @f$ a!\,b!\,c!/(a+b+c+d)! @f$.
   * - Wedge (Triangle @f$ \times [0,1] @f$): simplex moment times
   *   @f$ 1/(c+1) @f$.
   * - Pyramid @f$ \{0\le z\le 1,\ 0\le x,y\le 1-z\} @f$, the domain swept by
   *   the collapse @f$ (x,y,z)=((1-z)u,(1-z)v,z) @f$ used in
   *   GaussLegendre::buildPyramid:
   *   @f$ c!\,(a+b+2)!\,/\,\big[(a+b+c+3)!\,(a+1)(b+1)\big] @f$.
   */
  inline Real exactMoment(Geometry::Polytope::Type g, size_t a, size_t b, size_t c)
  {
    using Type = Geometry::Polytope::Type;
    switch (g)
    {
      case Type::Point:
        return (a == 0 && b == 0 && c == 0) ? Real(1) : Real(0);
      case Type::Segment:
        return Real(1) / static_cast<Real>(a + 1);
      case Type::Quadrilateral:
        return Real(1) / static_cast<Real>((a + 1) * (b + 1));
      case Type::Hexahedron:
        return Real(1) / static_cast<Real>((a + 1) * (b + 1) * (c + 1));
      case Type::Triangle:
        return factorial(a) * factorial(b) / factorial(a + b + 2);
      case Type::Tetrahedron:
        return factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3);
      case Type::Wedge:
        return (factorial(a) * factorial(b) / factorial(a + b + 2)) /
          static_cast<Real>(c + 1);
      case Type::Pyramid:
        return factorial(c) * factorial(a + b + 2) /
          (factorial(a + b + c + 3) * static_cast<Real>((a + 1) * (b + 1)));
    }
    return std::numeric_limits<Real>::quiet_NaN();
  }

  /// @brief Measure of the reference polytope, i.e. @f$ \int_K 1 \, dK @f$.
  inline Real referenceMeasure(Geometry::Polytope::Type g)
  {
    return exactMoment(g, 0, 0, 0);
  }

  /// @brief Spatial dimension of the reference polytope.
  /// @see Rodin::Geometry::Polytope::Traits::getDimension
  inline size_t dimensionOf(Geometry::Polytope::Type g)
  {
    return Geometry::Polytope::Traits(g).getDimension();
  }

  /// @brief Outcome of an exactness sweep over all monomials of a given degree.
  struct ExactnessReport
  {
      Real worstAbsoluteError = 0;   ///< Largest @f$ |Q(m) - I(m)| @f$ observed.
      Real worstRelativeError = 0;   ///< Largest error relative to @f$ |I(m)| @f$.
      size_t worstExponents[3] = {0, 0, 0}; ///< Monomial attaining the worst error.
      size_t monomialsTested = 0;    ///< Number of monomials in the sweep.
  };

  /**
   * @brief Integrates every monomial of total degree @f$ \le p @f$ and reports
   * the worst departure from the exact moment.
   *
   * This is the table-integrity instrument: a rule whose coefficients are
   * corrupted, mistyped, or solved to the wrong branch fails here, whereas
   * inspecting the numbers by eye does not.
   */
  template <class Rule>
  ExactnessReport exactnessSweep(const Rule& qf, Geometry::Polytope::Type g, size_t p)
  {
    ExactnessReport report;
    const size_t d = dimensionOf(g);
    const size_t amax = (d >= 1) ? p : 0;
    for (size_t a = 0; a <= amax; ++a)
    {
      const size_t bmax = (d >= 2) ? p - a : 0;
      for (size_t b = 0; b + a <= p && b <= bmax; ++b)
      {
        const size_t cmax = (d >= 3) ? p - a - b : 0;
        for (size_t c = 0; c + b + a <= p && c <= cmax; ++c)
        {
          Real numeric = 0;
          for (size_t i = 0; i < qf.getSize(); ++i)
          {
            const auto& x = qf.getPoint(i);
            Real m = 1;
            if (d >= 1)
              m *= std::pow(x[0], static_cast<Real>(a));
            if (d >= 2)
              m *= std::pow(x[1], static_cast<Real>(b));
            if (d >= 3)
              m *= std::pow(x[2], static_cast<Real>(c));
            numeric += qf.getWeight(i) * m;
          }
          const Real exact = exactMoment(g, a, b, c);
          const Real abserr = std::abs(numeric - exact);
          const Real relerr = abserr / std::max(std::abs(exact), Real(1e-300));
          ++report.monomialsTested;
          if (relerr > report.worstRelativeError)
          {
            report.worstRelativeError = relerr;
            report.worstAbsoluteError = abserr;
            report.worstExponents[0] = a;
            report.worstExponents[1] = b;
            report.worstExponents[2] = c;
          }
        }
      }
    }
    return report;
  }

  /// @brief Sum of the weights.
  template <class Rule>
  Real weightSum(const Rule& qf)
  {
    Real s = 0;
    for (size_t i = 0; i < qf.getSize(); ++i)
      s += qf.getWeight(i);
    return s;
  }

  /**
   * @brief Amplification factor @f$ \sum_i |w_i| / \sum_i w_i @f$.
   *
   * Equal to one exactly when every weight is non-negative. It bounds both the
   * roundoff amplification of the rule and the constant in the Peano-type
   * error estimate @f$ |E(f)| \le (|K| + \sum|w_i|)\,E_p(f) @f$, so it is the
   * single number separating a positive rule from a signed one.
   */
  template <class Rule>
  Real weightAmplification(const Rule& qf)
  {
    Real s = 0, as = 0;
    for (size_t i = 0; i < qf.getSize(); ++i)
    {
      s += qf.getWeight(i);
      as += std::abs(qf.getWeight(i));
    }
    return as / s;
  }

  /// @brief Whether every weight is strictly positive.
  template <class Rule>
  bool allWeightsPositive(const Rule& qf)
  {
    for (size_t i = 0; i < qf.getSize(); ++i)
      if (!(qf.getWeight(i) > Real(0)))
        return false;
    return true;
  }

  /**
   * @brief Whether @p x lies in the reference polytope @p g, up to @p tol.
   *
   * Uses the half-space description @f$ Ax \le b @f$ published by
   * Rodin::Geometry::Polytope::Traits, rather than a restatement of it: the
   * reference elements are defined there and nowhere else, so a test that
   * hard-codes them would pass against its own assumptions instead of against
   * the library.
   */
  inline bool isInside(
    Geometry::Polytope::Type g, const Math::SpatialVector<Real>& x, Real tol)
  {
    const Geometry::Polytope::Traits traits(g);
    const size_t d = traits.getDimension();
    if (d == 0)
      return true;
    const auto& hs = traits.getHalfSpace();
    Math::Vector<Real> p(static_cast<Eigen::Index>(d));
    for (size_t i = 0; i < d; ++i)
      p(static_cast<Eigen::Index>(i)) = x[i];
    const Math::Vector<Real> r = hs.matrix * p - hs.vector;
    for (Eigen::Index i = 0; i < r.size(); ++i)
      if (r(i) > tol)
        return false;
    return true;
  }

  /// @brief Whether every point lies in the reference polytope.
  template <class Rule>
  bool allPointsInside(const Rule& qf, Geometry::Polytope::Type g, Real tol = 1e-13)
  {
    for (size_t i = 0; i < qf.getSize(); ++i)
      if (!isInside(g, qf.getPoint(i), tol))
        return false;
    return true;
  }

  /// @brief An affine map @f$ x \mapsto A x + b @f$ of the reference element.
  struct AffineMap
  {
      Math::Matrix<Real> matrix; ///< The linear part @f$ A @f$.
      Math::Vector<Real> vector; ///< The translation @f$ b @f$.
  };

  /**
   * @brief The affine symmetry group of a reference polytope.
   *
   * The group is discovered rather than tabulated: every permutation of the
   * vertices is fitted with an affine map in the least-squares sense, and the
   * fit is kept only when it reproduces every vertex exactly. This admits
   * precisely the permutations realised by an affine self-map, so the same
   * code yields the six maps of the triangle, the twenty-four of the
   * tetrahedron, the forty-eight of the hexahedron, and the smaller groups of
   * the wedge and the pyramid, whose vertices are not interchangeable.
   */
  inline std::vector<AffineMap> symmetryGroup(
    Geometry::Polytope::Type g, Real tol = 1e-12)
  {
    const Geometry::Polytope::Traits traits(g);
    const auto d = static_cast<Eigen::Index>(traits.getDimension());
    const auto n = static_cast<Eigen::Index>(traits.getVertexCount());
    if (d == 0)
      return {};
    // Vertices in homogeneous form, so that one solve yields A and b together.
    Math::Matrix<Real> source(n, d + 1);
    for (Eigen::Index i = 0; i < n; ++i)
    {
      const auto& v = traits.getVertex(static_cast<size_t>(i));
      for (Eigen::Index k = 0; k < d; ++k)
        source(i, k) = v(k);
      source(i, d) = 1;
    }
    std::vector<size_t> permutation(static_cast<size_t>(n));
    std::iota(permutation.begin(), permutation.end(), size_t{0});
    const Eigen::FullPivHouseholderQR<Math::Matrix<Real>> qr(source);
    std::vector<AffineMap> group;
    do
    {
      Math::Matrix<Real> target(n, d);
      for (Eigen::Index i = 0; i < n; ++i)
      {
        const auto& v = traits.getVertex(permutation[static_cast<size_t>(i)]);
        for (Eigen::Index k = 0; k < d; ++k)
          target(i, k) = v(k);
      }
      const Math::Matrix<Real> fit = qr.solve(target);
      if ((source * fit - target).cwiseAbs().maxCoeff() > tol)
        continue;
      AffineMap map;
      map.matrix = fit.topRows(d).transpose();
      map.vector = fit.row(d).transpose();
      group.push_back(std::move(map));
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return group;
  }

  /**
   * @brief Whether a rule is carried onto itself by an affine map.
   *
   * Each image point must coincide with a node carrying the same weight, and
   * each node must be used once, so the nodes and weights together form an
   * invariant multiset rather than merely an invariant point set.
   */
  template <class Rule>
  bool isInvariantUnder(const Rule& qf, const AffineMap& map, Real tol = 1e-12)
  {
    const auto d = map.matrix.rows();
    const size_t size = qf.getSize();
    std::vector<bool> matched(size, false);
    for (size_t i = 0; i < size; ++i)
    {
      const auto& x = qf.getPoint(i);
      Math::Vector<Real> p(d);
      for (Eigen::Index k = 0; k < d; ++k)
        p(k) = x[static_cast<size_t>(k)];
      const Math::Vector<Real> image = map.matrix * p + map.vector;
      bool found = false;
      for (size_t j = 0; j < size && !found; ++j)
      {
        if (matched[j])
          continue;
        const auto& y = qf.getPoint(j);
        Real distance = 0;
        for (Eigen::Index k = 0; k < d; ++k)
          distance = std::max(distance, std::abs(image(k) - y[static_cast<size_t>(k)]));
        if (distance > tol)
          continue;
        if (std::abs(qf.getWeight(i) - qf.getWeight(j)) > tol)
          continue;
        matched[j] = true;
        found = true;
      }
      if (!found)
        return false;
    }
    return true;
  }

  /// @brief Whether a rule is invariant under the whole symmetry group.
  template <class Rule>
  bool isFullySymmetric(const Rule& qf, Geometry::Polytope::Type g, Real tol = 1e-12)
  {
    for (const auto& map : symmetryGroup(g))
      if (!isInvariantUnder(qf, map, tol))
        return false;
    return true;
  }
}

#endif
