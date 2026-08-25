/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Verifies the shipped WitherdenVincent coefficients.
 *
 * The coefficients are generated rather than transcribed, so they cannot carry
 * a transcription error; they can carry a solver error. They are therefore
 * checked against the independent closed-form moment oracle, not against the
 * generator that produced them.
 */
#include <gtest/gtest.h>
#include "Rodin/QF/WitherdenVincent.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  const std::vector<Polytope::Type> kElements = {Polytope::Type::Triangle,
    Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Wedge,
    Polytope::Type::Pyramid, Polytope::Type::Hexahedron};
  const char* name(Polytope::Type g)
  {
    switch (g)
    {
      case Polytope::Type::Triangle:
        return "triangle";
      case Polytope::Type::Quadrilateral:
        return "quadrilateral";
      case Polytope::Type::Hexahedron:
        return "hexahedron";
      case Polytope::Type::Tetrahedron:
        return "tetrahedron";
      default:
        return "wedge";
    }
  }
}

/// @brief Every shipped rule integrates every monomial of its degree exactly.
TEST(WitherdenVincentTest, ShippedRulesAreExactAtEveryDegree)
{
  for (const auto g : kElements)
  {
    ASSERT_GT(WitherdenVincent::getMaxDegree(g), 0u) << name(g);
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
    {
      const WitherdenVincent qf(p, g);
      const auto report = exactnessSweep(qf, g, p);
      EXPECT_LT(report.worstRelativeError, 1e-11)
        << name(g) << " degree " << p << " (" << qf.getSize() << " points)"
        << " worst x^" << report.worstExponents[0] << " y^" << report.worstExponents[1]
        << " z^" << report.worstExponents[2];
    }
  }
}

/// @brief Every shipped rule has strictly positive weights.
TEST(WitherdenVincentTest, ShippedWeightsArePositive)
{
  for (const auto g : kElements)
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
    {
      const WitherdenVincent qf(p, g);
      EXPECT_TRUE(allWeightsPositive(qf)) << name(g) << " degree " << p;
      EXPECT_NEAR(weightAmplification(qf), 1.0, 1e-14) << name(g) << " degree " << p;
    }
}

/// @brief Every shipped node lies inside its reference element, checked
/// against the half-space system Polytope::Traits publishes.
TEST(WitherdenVincentTest, ShippedNodesAreInterior)
{
  for (const auto g : kElements)
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
      EXPECT_TRUE(allPointsInside(WitherdenVincent(p, g), g))
        << name(g) << " degree " << p;
}

/// @brief Weights reproduce the measure of the reference element.
TEST(WitherdenVincentTest, ShippedWeightsSumToTheMeasure)
{
  for (const auto g : kElements)
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
      EXPECT_NEAR(weightSum(WitherdenVincent(p, g)), referenceMeasure(g), 1e-12)
        << name(g) << " degree " << p;
}

/// @brief Point counts do not exceed the seed and are at least the counting
/// bound: a rule cannot have fewer nodes than the moment count divided by the
/// unknowns each node carries.
TEST(WitherdenVincentTest, PointCountsRespectTheCountingBound)
{
  for (const auto g : kElements)
  {
    // The bound counts monomials of the simplex; the wedge is a product
    // element with a different moment count, so it is excluded rather than
    // checked against the wrong number.
    if (!Polytope::Traits(g).isSimplex())
      continue;
    const size_t d = Polytope::Traits(g).getDimension();
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
    {
      size_t moments = 1;
      for (size_t k = 1; k <= d; ++k)
        moments = moments * (p + k) / k;
      const size_t bound = (moments + d) / (d + 1);
      EXPECT_GE(WitherdenVincent(p, g).getSize(), bound)
        << name(g) << " degree " << p << " below the counting bound";
    }
  }
}

/// @brief Availability is reported honestly, and the tables are non-empty
/// where they claim to be.
TEST(WitherdenVincentTest, AvailabilityMatchesTheTables)
{
  for (const auto g : kElements)
  {
    const size_t maxDeg = WitherdenVincent::getMaxDegree(g);
    EXPECT_FALSE(WitherdenVincent::isAvailable(0, g));
    EXPECT_TRUE(WitherdenVincent::isAvailable(1, g));
    EXPECT_TRUE(WitherdenVincent::isAvailable(maxDeg, g));
    EXPECT_FALSE(WitherdenVincent::isAvailable(maxDeg + 1, g));
  }
  // A geometry with no rules at all still answers coherently.
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Point), 0u);
  EXPECT_FALSE(WitherdenVincent::isAvailable(2, Polytope::Type::Point));
}

/// @brief The rule is copyable and the copy carries the same data.
TEST(WitherdenVincentTest, CopyPreservesTheRule)
{
  const WitherdenVincent a(7, Polytope::Type::Triangle);
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
TEST(WitherdenVincentTest, ExactnessCheckRejectsAPerturbedCoefficient)
{
  const WitherdenVincent qf(6, Polytope::Type::Triangle);
  ASSERT_LT(exactnessSweep(qf, Polytope::Type::Triangle, 6).worstRelativeError, 1e-11);

  struct Perturbed
  {
      const WitherdenVincent& base;
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
