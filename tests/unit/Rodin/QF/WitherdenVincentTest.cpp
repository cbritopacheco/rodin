/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Verifies the shipped WitherdenVincent coefficients.
 *
 * The tabulated coefficients are checked against an independent closed-form
 * moment oracle and against geometric and positivity invariants.
 */
#include <gtest/gtest.h>
#include "PublishedCounts.h"
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
      case Polytope::Type::Wedge:
        return "wedge";
      case Polytope::Type::Pyramid:
        return "pyramid";
      default:
        return "unknown";
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

/**
 * @brief The discovered symmetry groups have their textbook orders.
 *
 * This guards the symmetry test against passing vacuously: were the search to
 * return only the identity, every rule would be trivially invariant. The
 * orders are those of the vertex symmetries of each reference element.
 */
TEST(WitherdenVincentTest, SymmetryGroupsHaveTheExpectedOrder)
{
  const std::vector<std::pair<Polytope::Type, size_t>> expected = {
    {Polytope::Type::Triangle, 6},   // S3
    {Polytope::Type::Quadrilateral, 8},   // D4
    {Polytope::Type::Tetrahedron, 24},  // S4
    {Polytope::Type::Wedge, 12},  // D3 x C2
    {Polytope::Type::Pyramid, 8},   // C4v
    {Polytope::Type::Hexahedron, 48}};  // the full octahedral group
  for (const auto& [g, order] : expected)
    EXPECT_EQ(symmetryGroup(g).size(), order) << name(g);
}

/**
 * @brief Every shipped rule is invariant under the element's symmetry group.
 *
 * Full symmetry is what distinguishes this family from
 * @ref Rodin::QF::XiaoGimbutas, and it is a property of the nodes and weights
 * themselves, not of the integrals they produce: an asymmetric rule can still
 * integrate exactly. The group is discovered from the reference vertices, and
 * hoisted out of the degree loop because discovering it is the expensive part.
 */
TEST(WitherdenVincentTest, ShippedRulesAreFullySymmetric)
{
  for (const auto g : kElements)
  {
    const auto group = symmetryGroup(g);
    ASSERT_GT(group.size(), 1u) << name(g) << " has no symmetries to check";
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
    {
      const WitherdenVincent qf(p, g);
      for (size_t s = 0; s < group.size(); ++s)
        EXPECT_TRUE(isInvariantUnder(qf, group[s]))
          << name(g) << " degree " << p << " is not invariant under symmetry " << s;
    }
  }
}

/**
 * @brief The symmetry check rejects a rule whose nodes have been displaced.
 *
 * Moving a single node off its orbit leaves the rule's size and its element
 * untouched, so this pins the check to the node positions.
 */
TEST(WitherdenVincentTest, SymmetryCheckRejectsADisplacedNode)
{
  struct Displaced
  {
      size_t getSize() const
      {
        return m_qf.getSize();
      }

      Math::SpatialVector<Real> getPoint(size_t i) const
      {
        auto x = m_qf.getPoint(i);
        if (i == 0)
          x[0] += 1e-3;
        return x;
      }

      Real getWeight(size_t i) const
      {
        return m_qf.getWeight(i);
      }

      const WitherdenVincent& m_qf;
  };

  const auto g = Polytope::Type::Triangle;
  const WitherdenVincent qf(4, g);
  ASSERT_TRUE(isFullySymmetric(qf, g));
  EXPECT_FALSE(isFullySymmetric(Displaced{qf}, g));
}

/// @brief Weights reproduce the measure of the reference element.
TEST(WitherdenVincentTest, ShippedWeightsSumToTheMeasure)
{
  for (const auto g : kElements)
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
      EXPECT_NEAR(weightSum(WitherdenVincent(p, g)), referenceMeasure(g), 1e-12)
        << name(g) << " degree " << p;
}

/// @brief Availability is reported honestly, and the tables are non-empty
/// where they claim to be.
TEST(WitherdenVincentTest, AvailabilityMatchesTheTables)
{
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Triangle), 20u);
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Quadrilateral), 20u);
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Tetrahedron), 10u);
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Wedge), 10u);
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Pyramid), 10u);
  EXPECT_EQ(WitherdenVincent::getMaxDegree(Polytope::Type::Hexahedron), 10u);
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

/**
 * @brief Every tabulated rule meets the published point count.
 *
 * This is the claim the tables exist to make, so it is asserted rather than
 * left to a changelog: a regression that quietly ships a larger rule would
 * otherwise pass every other test here, since a larger rule is still exact.
 *
 * Exactness and geometric validity are checked independently above.
 */
TEST(WitherdenVincentTest, TabulatedCountsMeetThePublishedOnes)
{
  for (const auto g : kElements)
  {
    for (size_t degree = 1; degree <= WitherdenVincent::getMaxDegree(g); ++degree)
    {
      const WitherdenVincent rule(degree, g);
      const size_t published = Rodin::Tests::QF::publishedCount(
        Rodin::Tests::QF::witherdenVincentCounts(g), degree);
      if (published == 0)
        continue;
      EXPECT_EQ(rule.getSize(), published)
        << name(g) << ", strength " << degree
        << ": point count differs from the published rule";
    }
  }
}
