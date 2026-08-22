/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <cmath>
#include "Rodin/QF/GaussLegendre.h"
using namespace Rodin;
using namespace Rodin::QF;

namespace
{
  /// @brief int_0^1 z^k (1-z)^a dz = k! a! / (k+a+1)!
  Real exactJacobiMoment(size_t k, size_t a)
  {
    const auto fact = [](size_t m) { Real r = 1; for (size_t i = 2; i <= m; ++i) r *= i; return r; };
    return fact(k) * fact(a) / fact(k + a + 1);
  }
}

/// @brief The Jacobi rule integrates z^k against (1-z)^alpha exactly up to
/// degree 2n-1. Verified against the closed-form Beta integral rather than
/// against another rule.
TEST(GaussJacobiTest, IsExactAgainstTheWeight)
{
  for (size_t alpha = 0; alpha <= 3; ++alpha)
    for (size_t n = 1; n <= 12; ++n)
    {
      std::vector<Real> x, w;
      GaussLegendre::gj1dUnit(n, alpha, x, w);
      ASSERT_EQ(x.size(), n);
      for (size_t k = 0; k <= 2 * n - 1; ++k)
      {
        Real s = 0;
        for (size_t i = 0; i < n; ++i)
          s += w[i] * std::pow(x[i], static_cast<Real>(k));
        const Real e = exactJacobiMoment(k, alpha);
        EXPECT_NEAR(s, e, 1e-13 * std::max(Real(1), std::abs(e)))
          << "alpha " << alpha << " n " << n << " k " << k;
      }
    }
}

/// @brief Weights are positive and nodes lie strictly inside (0,1).
TEST(GaussJacobiTest, WeightsPositiveAndNodesInterior)
{
  for (size_t alpha = 0; alpha <= 3; ++alpha)
    for (size_t n = 1; n <= 12; ++n)
    {
      std::vector<Real> x, w;
      GaussLegendre::gj1dUnit(n, alpha, x, w);
      for (size_t i = 0; i < n; ++i)
      {
        EXPECT_GT(w[i], 0) << "alpha " << alpha << " n " << n;
        EXPECT_GT(x[i], 0);
        EXPECT_LT(x[i], 1);
      }
    }
}

/// @brief With alpha = 0 the rule must reduce to Gauss-Legendre on [0,1],
/// whose nodes on a segment the dispatch already exposes.
TEST(GaussJacobiTest, AlphaZeroReducesToGaussLegendre)
{
  for (size_t n = 1; n <= 10; ++n)
  {
    std::vector<Real> xj, wj;
    GaussLegendre::gj1dUnit(n, 0, xj, wj);
    const GaussLegendre gl(Geometry::Polytope::Type::Segment, n);
    ASSERT_EQ(gl.getSize(), n);
    std::vector<Real> xl;
    for (size_t i = 0; i < n; ++i)
      xl.push_back(gl.getPoint(i)[0]);
    std::sort(xj.begin(), xj.end());
    std::sort(xl.begin(), xl.end());
    for (size_t i = 0; i < n; ++i)
      EXPECT_NEAR(xj[i], xl[i], 1e-13) << "n " << n << " node " << i;
  }
}

/// @brief The test itself must be able to fail: a deliberately wrong weight
/// exponent is rejected.
TEST(GaussJacobiTest, DetectsTheWrongWeight)
{
  std::vector<Real> x, w;
  GaussLegendre::gj1dUnit(4, 2, x, w);
  Real s = 0;
  for (size_t i = 0; i < x.size(); ++i)
    s += w[i];
  EXPECT_NEAR(s, exactJacobiMoment(0, 2), 1e-14);
  EXPECT_GT(std::abs(s - exactJacobiMoment(0, 1)), 1e-3)
    << "alpha=2 rule must not integrate the alpha=1 weight";
}
