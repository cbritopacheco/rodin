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
