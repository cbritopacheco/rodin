/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Verifies the shipped XiaoGimbutas coefficients.
 *
 * The tabulated coefficients are checked against an independent closed-form
 * moment oracle and against geometric and positivity invariants.
 */
#include <gtest/gtest.h>
#include "PublishedCounts.h"
#include "Rodin/QF/XiaoGimbutas.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  const std::vector<Polytope::Type> kSimplices = {
    Polytope::Type::Triangle, Polytope::Type::Tetrahedron};
  const char* name(Polytope::Type g)
  {
    return g == Polytope::Type::Triangle ? "triangle" : "tetrahedron";
  }
}

/// @brief Every shipped rule integrates every monomial of its degree exactly.
TEST(XiaoGimbutasTest, ShippedRulesAreExactAtEveryDegree)
{
  for (const auto g : kSimplices)
  {
    ASSERT_GT(XiaoGimbutas::getMaxDegree(g), 0u) << name(g);
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
    {
      const XiaoGimbutas qf(p, g);
      const auto report = exactnessSweep(qf, g, p);
      EXPECT_LT(report.worstRelativeError, 1e-11)
        << name(g) << " degree " << p << " (" << qf.getSize() << " points)"
        << " worst x^" << report.worstExponents[0] << " y^" << report.worstExponents[1]
        << " z^" << report.worstExponents[2];
    }
  }
}

/// @brief Every shipped rule has strictly positive weights.
TEST(XiaoGimbutasTest, ShippedWeightsArePositive)
{
  for (const auto g : kSimplices)
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
    {
      const XiaoGimbutas qf(p, g);
      EXPECT_TRUE(allWeightsPositive(qf)) << name(g) << " degree " << p;
      EXPECT_NEAR(weightAmplification(qf), 1.0, 1e-14) << name(g) << " degree " << p;
    }
}

/// @brief Every shipped node lies inside its reference element, checked
/// against the half-space system Polytope::Traits publishes.
TEST(XiaoGimbutasTest, ShippedNodesAreInterior)
{
  for (const auto g : kSimplices)
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
      EXPECT_TRUE(allPointsInside(XiaoGimbutas(p, g), g)) << name(g) << " degree " << p;
}

/**
 * @brief The shipped rules are not fully symmetric, and that is the point.
 *
 * Node elimination buys a smaller rule by giving up invariance under the
 * element's symmetry group, which is the trade this family makes against
 * @ref Rodin::QF::WitherdenVincent. Asserting the asymmetry keeps the two
 * families honestly distinguished, and confirms that the symmetry check used
 * on the symmetric family can in fact fail.
 */
TEST(XiaoGimbutasTest, ShippedRulesAreNotFullySymmetric)
{
  for (const auto g : kSimplices)
  {
    size_t asymmetric = 0;
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
      if (!isFullySymmetric(XiaoGimbutas(p, g), g))
        ++asymmetric;
    EXPECT_GT(asymmetric, 0u)
      << name(g) << ": every shipped rule is symmetric, so nothing was gained"
      << " over the symmetric family";
  }
}

/// @brief Weights reproduce the measure of the reference element.
TEST(XiaoGimbutasTest, ShippedWeightsSumToTheMeasure)
{
  for (const auto g : kSimplices)
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
      EXPECT_NEAR(weightSum(XiaoGimbutas(p, g)), referenceMeasure(g), 1e-12)
        << name(g) << " degree " << p;
}

/// @brief Availability is reported honestly, and the tables are non-empty
/// where they claim to be.
TEST(XiaoGimbutasTest, AvailabilityMatchesTheTables)
{
  EXPECT_EQ(XiaoGimbutas::getMaxDegree(Polytope::Type::Triangle), 50u);
  EXPECT_EQ(XiaoGimbutas::getMaxDegree(Polytope::Type::Tetrahedron), 15u);
  for (const auto g : kSimplices)
  {
    const size_t maxDeg = XiaoGimbutas::getMaxDegree(g);
    EXPECT_FALSE(XiaoGimbutas::isAvailable(0, g));
    EXPECT_TRUE(XiaoGimbutas::isAvailable(1, g));
    EXPECT_TRUE(XiaoGimbutas::isAvailable(maxDeg, g));
    EXPECT_FALSE(XiaoGimbutas::isAvailable(maxDeg + 1, g));
  }
  EXPECT_EQ(XiaoGimbutas::getMaxDegree(Polytope::Type::Hexahedron), 0u);
  EXPECT_FALSE(XiaoGimbutas::isAvailable(2, Polytope::Type::Hexahedron));
}

/// @brief The rule is copyable and the copy carries the same data.
TEST(XiaoGimbutasTest, CopyPreservesTheRule)
{
  const XiaoGimbutas a(7, Polytope::Type::Triangle);
  const std::unique_ptr<QuadratureFormulaBase> b(a.copy());
  ASSERT_EQ(b->getSize(), a.getSize());
  for (size_t i = 0; i < a.getSize(); ++i)
  {
    EXPECT_EQ(b->getWeight(i), a.getWeight(i));
    for (Eigen::Index k = 0; k < a.getPoint(i).size(); ++k)
      EXPECT_EQ(b->getPoint(i)[k], a.getPoint(i)[k]);
  }
}

/// @brief The exactness check applied here can fail: perturbing one shipped
/// weight must be rejected.
TEST(XiaoGimbutasTest, ExactnessCheckRejectsAPerturbedCoefficient)
{
  const XiaoGimbutas qf(6, Polytope::Type::Triangle);
  ASSERT_LT(exactnessSweep(qf, Polytope::Type::Triangle, 6).worstRelativeError, 1e-11);

  struct Perturbed
  {
      const XiaoGimbutas& base;
      size_t getSize() const
      {
        return base.getSize();
      }
      Real getWeight(size_t i) const
      {
        return base.getWeight(i) * (i == 0 ? 1.0 + 1e-9 : 1.0);
      }
      const Math::SpatialVector<Real>& getPoint(size_t i) const
      {
        return base.getPoint(i);
      }
  } bad{qf};

  EXPECT_GT(exactnessSweep(bad, Polytope::Type::Triangle, 6).worstRelativeError, 1e-11);
}

/**
 * @brief Every tabulated rule has the published point count.
 *
 * These rules are measured against the Xiao--Gimbutas counts rather than the
 * symmetric Witherden--Vincent counts. The two families publish different
 * numbers for the same element and cannot be compared interchangeably.
 */
TEST(XiaoGimbutasTest, TabulatedCountsMeetThePublishedOnes)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    for (size_t degree = 1; degree <= XiaoGimbutas::getMaxDegree(g); ++degree)
    {
      const XiaoGimbutas rule(degree, g);
      const size_t published =
        Rodin::Tests::QF::publishedCount(Rodin::Tests::QF::xiaoGimbutasCounts(g), degree);
      if (published == 0)
        continue;
      EXPECT_EQ(rule.getSize(), published)
        << name(g) << ", strength " << degree
        << ": point count differs from the published rule";
    }
  }
}
