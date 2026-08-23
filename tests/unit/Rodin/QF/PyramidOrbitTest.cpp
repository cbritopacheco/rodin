/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Checks the D4 pyramid orbits before any rule is solved with them.
 *
 * The shear onto the corner-apex reference element is the kind of transform
 * whose errors are silent: a wrong one still produces plausible points inside
 * a plausible element. It is verified here on its own.
 */
#include <gtest/gtest.h>
#include <random>
#include "Rodin/QF/SymmetricOrbit.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

using PC = SymmetricOrbit::PyramidClass;

/// @brief Cardinalities match the classes they name.
TEST(PyramidOrbitTest, ClassCardinalities)
{
  EXPECT_EQ(SymmetricOrbit::pyramidClassSize(PC::Centre), 1u);
  EXPECT_EQ(SymmetricOrbit::pyramidClassSize(PC::Axis), 4u);
  EXPECT_EQ(SymmetricOrbit::pyramidClassSize(PC::Diagonal), 4u);
  EXPECT_EQ(SymmetricOrbit::pyramidClassSize(PC::General), 8u);
  EXPECT_EQ(SymmetricOrbit::expandPyramid(PC::Centre, 0, 0, 0.3).size(), 1u);
  EXPECT_EQ(SymmetricOrbit::expandPyramid(PC::Axis, 0.4, 0, 0.3).size(), 4u);
  EXPECT_EQ(SymmetricOrbit::expandPyramid(PC::Diagonal, 0.4, 0, 0.3).size(), 4u);
  EXPECT_EQ(SymmetricOrbit::expandPyramid(PC::General, 0.5, 0.2, 0.3).size(), 8u);
}

/// @brief Every orbit point lands inside the reference pyramid, checked
/// against the half-space system Polytope::Traits publishes rather than
/// against my derivation of it.
TEST(PyramidOrbitTest, PointsLandInsideTheReferencePyramid)
{
  std::mt19937 rng(3u);
  std::uniform_real_distribution<Real> par(-0.95, 0.95), hgt(0.02, 0.95);
  for (const auto c : {PC::Centre, PC::Axis, PC::Diagonal, PC::General})
    for (int trial = 0; trial < 500; ++trial)
    {
      const Real a = par(rng), b = par(rng), z = hgt(rng);
      for (const auto& p : SymmetricOrbit::expandPyramid(c, a, b, z))
        ASSERT_TRUE(isInside(Polytope::Type::Pyramid, p, 1e-12))
          << "alpha " << a << " beta " << b << " z " << z;
    }
}

/// @brief The orbit is symmetric in centred coordinates: undoing the shear
/// must give a point set invariant under the dihedral group of the square.
/// This is what makes it an orbit rather than an arbitrary set of points.
TEST(PyramidOrbitTest, OrbitIsDihedralInCentredCoordinates)
{
  const Real z = 0.4, half = (1 - z) / 2;
  for (const auto c : {PC::Axis, PC::Diagonal, PC::General})
  {
    const auto pts = SymmetricOrbit::expandPyramid(c, 0.6, 0.25, z);
    std::vector<std::pair<Real, Real>> uv;
    for (const auto& p : pts)
      uv.emplace_back((p[0] - half) / half, (p[1] - half) / half);

    const auto contains = [&](Real u, Real v)
    {
      for (const auto& [a, b] : uv)
        if (std::abs(a - u) < 1e-12 && std::abs(b - v) < 1e-12)
          return true;
      return false;
    };
    for (const auto& [u, v] : uv)
    {
      EXPECT_TRUE(contains(-u, v)) << "not closed under u -> -u";
      EXPECT_TRUE(contains(u, -v)) << "not closed under v -> -v";
      EXPECT_TRUE(contains(v, u)) << "not closed under the diagonal swap";
    }
  }
}

/// @brief The shear has unit determinant, so it moves points without changing
/// any volume. Verified by mapping the reference pyramid's own vertices.
TEST(PyramidOrbitTest, ShearPreservesVolume)
{
  // At height z the cross-section is a square of side 1 - z in both the
  // centred and the reference description, so the swept volume agrees.
  for (const Real z : {0.0, 0.25, 0.5, 0.9})
  {
    const auto corners = SymmetricOrbit::expandPyramid(PC::Diagonal, 1.0, 0, z);
    Real minx = 1e300, maxx = -1e300;
    for (const auto& p : corners)
    {
      minx = std::min(minx, p[0]);
      maxx = std::max(maxx, p[0]);
    }
    EXPECT_NEAR(maxx - minx, 1 - z, 1e-14)
      << "cross-section width wrong at z = " << z;
    EXPECT_NEAR(minx, 0.0, 1e-14) << "cross-section not anchored at the corner";
  }
}

/// @brief Degenerate parameters do not escape the element: a zero shape
/// parameter collapses onto the axis, and the apex is the limit z -> 1.
TEST(PyramidOrbitTest, DegenerateParametersStayInside)
{
  for (const Real z : {0.0, 0.999})
    for (const auto c : {PC::Centre, PC::Axis, PC::Diagonal, PC::General})
      for (const auto& p : SymmetricOrbit::expandPyramid(c, 0, 0, z))
        EXPECT_TRUE(isInside(Polytope::Type::Pyramid, p, 1e-12));
}
