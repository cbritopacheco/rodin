/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Fixed-mesh degree convergence for vector linear elasticity. */

#include <string>

#include <gtest/gtest.h>

#include "LinearElasticity.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::Convergence;

namespace Rodin::Tests::Convergence::P::LinearElasticity
{
  class PConvergenceTest : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(PConvergenceTest, AnalyticDisplacementDecaysAtEveryDegree)
  {
    UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(2);
    const ::Rodin::Tests::Convergence::LinearElasticity::ManufacturedSolution
      data(grid.getDimension(), 1.5, 0.5);
    ErrorHistory history;
    history.append(1, ::Rodin::Tests::Convergence::LinearElasticity::solve<1>(mesh, data))
           .append(2, ::Rodin::Tests::Convergence::LinearElasticity::solve<2>(mesh, data))
           .append(3, ::Rodin::Tests::Convergence::LinearElasticity::solve<3>(mesh, data))
           .append(4, ::Rodin::Tests::Convergence::LinearElasticity::solve<4>(mesh, data));

    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      SCOPED_TRACE(::testing::Message()
        << "degrees " << i << " -> " << i + 1 << ": L2 "
        << coarse.getL2() << " -> " << fine.getL2()
        << ", H1 " << coarse.getH1Seminorm() << " -> "
        << fine.getH1Seminorm());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getExponentialRates(i);
      EXPECT_GT(rate.getL2(), 0.25);
      EXPECT_GT(rate.getH1Seminorm(), 0.25);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PConvergenceTest,
    ::testing::Values(
      Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    { return std::string(UniformGrid::getGeometryName(info.param)); });
}
