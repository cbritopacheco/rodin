/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Combined mesh/degree improvement for mixed Stokes fields. */

#include "../../StokesProblem.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::HP::Stokes
{
  class StokesHPTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(StokesHPTest, BothFieldsImproveAlongCombinedPath)
  {
    const UniformGrid grid(GetParam());
    const StokesData data(grid.getDimension(), StokesData::Field::Smooth);
    ErrorHistory velocity, pressure;
    const auto append = [&](Real h, const StokesErrors& errors) {
      velocity.append(h, errors.velocity);
      pressure.append(h, errors.pressure);
    };
    append(0.5, StokesProblem(grid.makeMesh(3), data).solve<2>());
    append(Real(1) / 3, StokesProblem(grid.makeMesh(4), data).solve<3>());
    append(0.25, StokesProblem(grid.makeMesh(5), data).solve<4>());
    for (bool isVelocity : {true, false})
    {
      const auto& history = isVelocity ? velocity : pressure;
      ASSERT_EQ(history.getSize(), 3u);
      for (size_t i = 1; i < history.getSize(); ++i)
      {
        const auto& coarse = history.getSample(i - 1).error;
        const auto& fine = history.getSample(i).error;
        SCOPED_TRACE(::testing::Message()
          << "velocity=" << isVelocity << " interval=" << i << " L2 " << coarse.getL2()
          << " -> " << fine.getL2() << " H1 " << coarse.getH1Seminorm() << " -> "
          << fine.getH1Seminorm());
        ASSERT_TRUE(coarse.isFinite());
        ASSERT_TRUE(fine.isFinite());
        ASSERT_GT(fine.getL2(), 0);
        ASSERT_GT(fine.getH1Seminorm(), 0);
        ASSERT_GT(coarse.getL2(), fine.getL2());
        ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
        const auto rate = history.getAlgebraicRates(i);
        SCOPED_TRACE(::testing::Message()
          << "effective rates " << rate.getL2() << ", " << rate.getH1Seminorm());
        EXPECT_GT(rate.getL2(), isVelocity ? 2.5 : 1.5);
        EXPECT_GT(rate.getH1Seminorm(), isVelocity ? 1.5 : 0.5);
      }
    }
  }

  TEST_P(StokesHPTest, P4P3ReproducesQuarticPatch)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const auto errors =
      StokesProblem(mesh, StokesData(mesh.getDimension(), StokesData::Field::Quartic))
        .solve<4>();
    EXPECT_LT(errors.velocity.getL2(), 1e-9);
    EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
    EXPECT_LT(errors.pressure.getL2(), 1e-9);
    EXPECT_LT(errors.pressure.getH1Seminorm(), 1e-9);
    EXPECT_LT(errors.divergence, 1e-9);
  }

  TEST_P(StokesHPTest, P4P3RejectsIncorrectViscosity)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const auto errors =
      StokesProblem(mesh, StokesData(mesh.getDimension(), StokesData::Field::Quadratic))
        .solve<4>(2);
    EXPECT_LT(errors.velocity.getL2(), 1e-9);
    EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
    EXPECT_GT(errors.pressure.getL2(), 0.1);
    EXPECT_GT(errors.pressure.getH1Seminorm(), 1);
  }

  TEST_P(StokesHPTest, P4QuadratureSensitivity)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const StokesData data(mesh.getDimension(), StokesData::Field::Smooth);
    const auto baseline = StokesProblem(mesh, data).solve<4>();
    const auto refined = StokesProblem(mesh, data, 18).solve<4>();
    for (const auto& pair : {std::pair{baseline.velocity, refined.velocity},
           std::pair{baseline.pressure, refined.pressure}})
    {
      ASSERT_GT(pair.first.getL2(), 0);
      ASSERT_GT(pair.first.getH1Seminorm(), 0);
      EXPECT_LT(std::abs(pair.second.getL2() / pair.first.getL2() - 1), 1e-6);
      EXPECT_LT(
        std::abs(pair.second.getH1Seminorm() / pair.first.getH1Seminorm() - 1), 1e-6);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, StokesHPTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
