/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Combined mesh and degree refinement for linear elasticity. */

#include <string>

#include <gtest/gtest.h>

#include "LinearElasticity.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::Convergence;

namespace Rodin::Tests::Convergence::HP::LinearElasticity
{
  class HPConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(HPConvergenceTest, AnalyticDisplacementImprovesAlongCombinedPath)
  {
    UniformGrid grid(GetParam());
    const ::Rodin::Tests::Convergence::LinearElasticity::ManufacturedSolution data(
      grid.getDimension(), 1.5, 0.5);
    ErrorHistory history;
    history
      .append(1,
        ::Rodin::Tests::Convergence::LinearElasticity::solve<1>(grid.makeMesh(2), data))
      .append(0.5,
        ::Rodin::Tests::Convergence::LinearElasticity::solve<2>(grid.makeMesh(3), data))
      .append(0.25,
        ::Rodin::Tests::Convergence::LinearElasticity::solve<3>(grid.makeMesh(5), data));

    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      SCOPED_TRACE(::testing::Message()
        << "interval " << i << ": L2 " << coarse.getL2() << " -> " << fine.getL2()
        << ", H1 " << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      EXPECT_GT(rate.getL2(), 1.9);
      EXPECT_GT(rate.getH1Seminorm(), 0.9);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, HPConvergenceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
