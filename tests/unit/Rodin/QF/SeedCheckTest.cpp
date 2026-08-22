/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"
using namespace Rodin; using namespace Rodin::QF;
using namespace Rodin::Geometry; using namespace Rodin::Tests::QF;

TEST(NodeEliminationTest, ProductSeedIsExactPositiveAndInterior)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
    for (size_t p = 1; p <= 12; ++p)
    {
      const auto seed = NodeElimination::productSeed(g, p);
      EXPECT_LT(exactnessSweep(seed, g, p).worstRelativeError, 1e-12)
        << "degree " << p << " n=" << seed.getSize();
      EXPECT_TRUE(allWeightsPositive(seed));
      EXPECT_TRUE(allPointsInside(seed, g));
    }
}

TEST(NodeEliminationTest, ReduceKeepsExactnessAndShrinks)
{
  for (size_t p = 1; p <= 8; ++p)
  {
    const auto seed = NodeElimination::productSeed(Polytope::Type::Triangle, p);
    const auto t0 = std::chrono::steady_clock::now();
    const auto out = NodeElimination::reduce(Polytope::Type::Triangle, p, seed);
    const double secs = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
    EXPECT_LE(out.getSize(), seed.getSize());
    EXPECT_LT(exactnessSweep(out, Polytope::Type::Triangle, p).worstRelativeError, 1e-11)
      << "degree " << p;
    EXPECT_TRUE(allWeightsPositive(out)) << "degree " << p;
    EXPECT_TRUE(allPointsInside(out, Polytope::Type::Triangle)) << "degree " << p;
    std::printf("deg %zu: seed %zu -> %zu points (%.1fs)\n",
      p, seed.getSize(), out.getSize(), secs);
  }
}
