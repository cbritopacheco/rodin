/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Verifies that the reference elements the quadrature work integrates
 * over are the ones Rodin defines.
 *
 * Three descriptions of each reference element are in play: the vertices and
 * the half-space system published by Geometry::Polytope::Traits, and the
 * closed-form moments the exactness oracle uses. A rule can be exact against
 * the wrong domain and every other test in the suite will pass, so the three
 * are pinned against each other here.
 */
#include <gtest/gtest.h>
#include <random>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  std::vector<Polytope::Type> elements()
  {
    std::vector<Polytope::Type> out;
    for (const auto g : Polytope::Types)
      if (g != Polytope::Type::Point)
        out.push_back(g);
    return out;
  }
}

/// @brief Every reference vertex satisfies the half-space system with zero
/// slack on at least one row: vertices lie on the boundary, not merely inside.
TEST(ReferenceElementTest, VerticesLieOnTheHalfSpaceBoundary)
{
  for (const auto g : elements())
  {
    const Polytope::Traits traits(g);
    const auto& hs = traits.getHalfSpace();
    const size_t d = traits.getDimension();
    for (size_t i = 0; i < traits.getVertexCount(); ++i)
    {
      Math::Vector<Real> v(static_cast<Eigen::Index>(d));
      for (size_t k = 0; k < d; ++k)
        v(static_cast<Eigen::Index>(k)) =
          traits.getVertex(i)[static_cast<Eigen::Index>(k)];
      const Math::Vector<Real> res = hs.matrix * v - hs.vector;
      EXPECT_LE(res.maxCoeff(), 1e-13) << "vertex " << i << " outside";
      EXPECT_GE(res.maxCoeff(), -1e-13) << "vertex " << i << " strictly inside";
    }
  }
}

/// @brief The half-space system encloses exactly the volume the moment oracle
/// reports. A missing or redundant constraint changes the region without
/// changing any vertex, and would otherwise go unnoticed.
TEST(ReferenceElementTest, HalfSpaceVolumeMatchesTheMomentOracle)
{
  std::mt19937 rng(20260101u);
  std::uniform_real_distribution<Real> uni(0, 1);
  for (const auto g : elements())
  {
    const Polytope::Traits traits(g);
    const size_t d = traits.getDimension();
    constexpr size_t samples = 400000;
    size_t hits = 0;
    for (size_t s = 0; s < samples; ++s)
    {
      Math::SpatialVector<Real> x;
      x.resize(static_cast<Eigen::Index>(d));
      for (size_t k = 0; k < d; ++k)
        x[static_cast<Eigen::Index>(k)] = uni(rng);
      if (isInside(g, x, 0))
        ++hits;
    }
    // The sampling box is the unit cube, whose volume is one, so the hit
    // fraction estimates the element's measure directly.
    const Real estimate = static_cast<Real>(hits) / static_cast<Real>(samples);
    const Real exact = referenceMeasure(g);
    EXPECT_NEAR(estimate, exact, 0.01)
      << "half-space region and moment oracle disagree on the measure";
  }
}

/// @brief The seed rule integrates over the same region: its points satisfy
/// the half-space system and its weights reproduce the oracle's measure.
TEST(ReferenceElementTest, ProductSeedAgreesWithTheReferenceElement)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    for (size_t p = 1; p <= 8; ++p)
    {
      const auto seed = NodeElimination::productSeed(g, p);
      EXPECT_TRUE(allPointsInside(seed, g)) << "degree " << p;
      EXPECT_NEAR(weightSum(seed), referenceMeasure(g), 1e-13) << "degree " << p;
    }
  }
}

/// @brief The admissibility predicate used inside the solver and the one used
/// by the tests accept the same points. They are written separately, and a
/// solver that policed a different region than the tests check would produce
/// rules that pass every test and are wrong.
TEST(ReferenceElementTest, SolverAndTestAdmissibilityAgree)
{
  std::mt19937 rng(7u);
  std::uniform_real_distribution<Real> uni(-0.2, 1.2);
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    const size_t d = Polytope::Traits(g).getDimension();
    size_t agreed = 0, inside = 0;
    for (size_t s = 0; s < 20000; ++s)
    {
      NodeElimination::Rule one;
      Math::SpatialVector<Real> x;
      x.resize(static_cast<Eigen::Index>(d));
      for (size_t k = 0; k < d; ++k)
        x[static_cast<Eigen::Index>(k)] = uni(rng);
      one.points.push_back(x);
      one.weights.push_back(1.0);

      const bool solver = NodeElimination::isAdmissible(g, one, 1e-9);
      const bool test = isInside(g, x, -1e-9);
      if (solver == test)
        ++agreed;
      if (test)
        ++inside;
    }
    EXPECT_EQ(agreed, 20000u) << "predicates disagree";
    EXPECT_GT(inside, 0u) << "sampling never landed inside; test is vacuous";
  }
}
