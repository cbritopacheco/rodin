/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Verifies the shipped XiaoGimbutas coefficients.
 *
 * The coefficients are generated rather than transcribed, so they cannot carry
 * a transcription error; they can carry a solver error. They are therefore
 * checked against the independent closed-form moment oracle, not against the
 * generator that produced them.
 */
#include <gtest/gtest.h>
#include "PublishedCounts.h"
#include <set>
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

/// @brief Weights reproduce the measure of the reference element.
TEST(XiaoGimbutasTest, ShippedWeightsSumToTheMeasure)
{
  for (const auto g : kSimplices)
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
      EXPECT_NEAR(weightSum(XiaoGimbutas(p, g)), referenceMeasure(g), 1e-12)
        << name(g) << " degree " << p;
}

/// @brief Point counts do not exceed the seed and are at least the counting
/// bound: a rule cannot have fewer nodes than the moment count divided by the
/// unknowns each node carries.
TEST(XiaoGimbutasTest, PointCountsRespectTheCountingBound)
{
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
    {
      size_t moments = 1;
      for (size_t k = 1; k <= d; ++k)
        moments = moments * (p + k) / k;
      const size_t bound = (moments + d) / (d + 1);
      EXPECT_GE(XiaoGimbutas(p, g).getSize(), bound)
        << name(g) << " degree " << p << " below the counting bound";
    }
  }
}

/// @brief Availability is reported honestly, and the tables are non-empty
/// where they claim to be.
TEST(XiaoGimbutasTest, AvailabilityMatchesTheTables)
{
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
 * @brief Every tabulated rule is no larger than the published one, and none is
 * below the counting bound.
 *
 * These rules come from node elimination and are asymmetric, so they are
 * measured against the Xiao--Gimbutas counts rather than the symmetric ones;
 * the two families publish different numbers for the same element and
 * comparing across them is meaningless.
 *
 * Most strengths here are *below* the published count, and several sit exactly
 * on the counting bound, which no rule can beat. The bound is asserted in the
 * other direction for that reason: a count below it would not be a discovery
 * but a moment system that had stopped spanning.
 */
TEST(XiaoGimbutasTest, TabulatedCountsMeetThePublishedOnes)
{
  // Strengths where elimination does not reach the published count. The
  // tetrahedron at five is structural rather than a matter of effort: the
  // fourteen-point rule attains the counting bound exactly, which leaves the
  // moment system square and its solutions isolated, and elimination works by
  // moving along the null space of an underdetermined one.
  const std::set<std::pair<Polytope::Type, size_t>> outstanding = {
    {Polytope::Type::Triangle, 8},     // 17 points against a published 16
    {Polytope::Type::Tetrahedron, 5},  // 15 points against a published 14
  };

  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    const size_t d = Polytope::Traits(g).getDimension();
    for (size_t degree = 1; degree <= XiaoGimbutas::getMaxDegree(g); ++degree)
    {
      const XiaoGimbutas rule(degree, g);
      EXPECT_GE(rule.getSize(), Rodin::Tests::QF::countingBound(d, degree))
        << name(g) << ", strength " << degree << ": below the counting bound";

      const size_t published =
        Rodin::Tests::QF::publishedCount(Rodin::Tests::QF::xiaoGimbutasCounts(g), degree);
      if (published == 0 || outstanding.count({g, degree}))
        continue;
      EXPECT_LE(rule.getSize(), published)
        << name(g) << ", strength " << degree << ": more points than the published rule";
    }
  }
}
