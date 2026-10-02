/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief H-convergence validation for the Taylor--Hood Stokes discretization.
 */

#include <cstdint>
#include <functional>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "../../StokesProblem.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::Stokes
{
  StokesErrors solve(
    const UniformGridHierarchy& hierarchy, size_t pointsPerAxis, const StokesData& data)
  {
    auto mesh = hierarchy.makeMesh(pointsPerAxis);
    return StokesProblem(mesh, data, 12).solve<2>();
  }

  void expectRates(const ErrorHistory& velocity, const ErrorHistory& pressure)
  {
    ASSERT_EQ(velocity.getSize(), pressure.getSize());
    ASSERT_GE(velocity.getSize(), 3);
    for (size_t i = 1; i < velocity.getSize(); ++i)
    {
      const auto& coarseVelocity = velocity.getSample(i - 1).error;
      const auto& fineVelocity = velocity.getSample(i).error;
      const auto& coarsePressure = pressure.getSample(i - 1).error;
      const auto& finePressure = pressure.getSample(i).error;
      ASSERT_TRUE(coarseVelocity.isFinite());
      ASSERT_TRUE(fineVelocity.isFinite());
      ASSERT_TRUE(coarsePressure.isFinite());
      ASSERT_TRUE(finePressure.isFinite());
      ASSERT_GT(coarseVelocity.getL2(), fineVelocity.getL2());
      ASSERT_GT(coarseVelocity.getH1Seminorm(), fineVelocity.getH1Seminorm());
      ASSERT_GT(coarsePressure.getL2(), finePressure.getL2());
      ASSERT_GT(coarsePressure.getH1Seminorm(), finePressure.getH1Seminorm());

      const auto velocityRates = velocity.getAlgebraicRates(i);
      const auto pressureRates = pressure.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "velocity L2 " << coarseVelocity.getL2() << " -> " << fineVelocity.getL2()
        << " (rate " << velocityRates.getL2() << "), velocity H1 "
        << coarseVelocity.getH1Seminorm() << " -> " << fineVelocity.getH1Seminorm()
        << " (rate " << velocityRates.getH1Seminorm() << "); pressure L2 "
        << coarsePressure.getL2() << " -> " << finePressure.getL2() << " (rate "
        << pressureRates.getL2() << "), pressure H1 " << coarsePressure.getH1Seminorm()
        << " -> " << finePressure.getH1Seminorm() << " (rate "
        << pressureRates.getH1Seminorm() << ')');
      EXPECT_GT(velocityRates.getL2(), 2.45);
      EXPECT_LT(velocityRates.getL2(), 3.55);
      EXPECT_GT(velocityRates.getH1Seminorm(), 1.55);
      EXPECT_LT(velocityRates.getH1Seminorm(), 2.45);
      EXPECT_GT(pressureRates.getL2(), 1.45);
      EXPECT_LT(pressureRates.getL2(), 2.55);
      EXPECT_GT(pressureRates.getH1Seminorm(), 0.55);
      EXPECT_LT(pressureRates.getH1Seminorm(), 1.45);
    }
  }

  class StokesHConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  /** @brief Verifies exact reproduction by the Taylor--Hood pair. */
  TEST_P(StokesHConvergenceTest, AffineDivergenceFreeSolutionIsExact)
  {
    const UniformGridHierarchy hierarchy(GetParam(), {3});
    const auto error = solve(
      hierarchy, 3, StokesData(hierarchy.getDimension(), StokesData::Field::Affine));
    EXPECT_LT(error.velocity.getL2(), 1e-10);
    EXPECT_LT(error.velocity.getH1Seminorm(), 1e-10);
    EXPECT_LT(error.pressure.getL2(), 1e-10);
    EXPECT_LT(error.pressure.getH1Seminorm(), 1e-10);
    EXPECT_LT(error.divergence, 1e-10);
  }

  /** @brief Verifies Taylor--Hood velocity and pressure h-rates. */
  TEST_P(StokesHConvergenceTest, PolynomialDivergenceFreeSolutionHasOptimalRates)
  {
    const UniformGridHierarchy hierarchy(GetParam(), {3, 5, 9});
    const StokesData data(hierarchy.getDimension(), StokesData::Field::Cubic);
    ErrorHistory velocity;
    ErrorHistory pressure;
    for (const size_t level : hierarchy.getLevels())
    {
      const auto error = solve(hierarchy, level, data);
      const Real h = hierarchy.getMeshSize(level);
      velocity.append(h, error.velocity);
      pressure.append(h, error.pressure);
    }
    expectRates(velocity, pressure);
  }

  std::string geometryName(const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries, StokesHConvergenceTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    geometryName);
}
