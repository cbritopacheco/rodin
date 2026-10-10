/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Semilinear Poisson combined-refinement verification. */

#include "../../NonlinearPoisson.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::HP::NonlinearPoisson
{
  class NonlinearPoissonHPTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(NonlinearPoissonHPTest, ErrorsImproveAlongCombinedPath)
  {
    const UniformGrid grid(GetParam());
    auto mesh1 = grid.makeMesh(2), mesh2 = grid.makeMesh(3), mesh3 = grid.makeMesh(5);
    ErrorHistory history;
    history.append(1, NonlinearPoissonProblem(mesh1).solve<1>());
    history.append(0.5, NonlinearPoissonProblem(mesh2).solve<2>());
    history.append(0.25, NonlinearPoissonProblem(mesh3).solve<3>());
    for (size_t i = 1; i < 3; ++i)
    {
      const auto& a = history.getSample(i - 1).error;
      const auto& b = history.getSample(i).error;
      ASSERT_TRUE(a.isFinite());
      ASSERT_TRUE(b.isFinite());
      ASSERT_GT(b.getL2(), 0);
      ASSERT_GT(b.getH1Seminorm(), 0);
      ASSERT_GT(a.getL2(), b.getL2());
      ASSERT_GT(a.getH1Seminorm(), b.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "interval=" << i << " L2 " << a.getL2() << " -> " << b.getL2() << " H1 "
        << a.getH1Seminorm() << " -> " << b.getH1Seminorm() << " effective rates "
        << rate.getL2() << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), 1.9);
      EXPECT_GT(rate.getH1Seminorm(), 0.9);
    }
  }

  TEST_P(NonlinearPoissonHPTest, TangentsAgreeWithResidualAndRejectWrongDerivative)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const NonlinearPoissonProblem problem(mesh);
    EXPECT_LT(problem.tangentDefect<1>(), 1e-6);
    EXPECT_LT(problem.tangentDefect<2>(), 1e-6);
    EXPECT_LT(problem.tangentDefect<3>(), 1e-6);
    EXPECT_GT(problem.tangentDefect<3>(true), 1e-3);
  }

  TEST_P(NonlinearPoissonHPTest, RejectsOmittedCubicReaction)
  {
    // Resolve the correct field below the same bounds used to reject the
    // incorrect physics, so discretization error cannot pass this control.
    auto mesh = UniformGrid(GetParam()).makeMesh(5);
    const NonlinearPoissonProblem problem(mesh, 16, 4);
    const auto baseline = problem.solve<3>();
    EXPECT_LT(baseline.getL2(), 0.05);
    EXPECT_LT(baseline.getH1Seminorm(), 0.2);
    const auto errors = problem.solve<3>(true);
    EXPECT_GT(errors.getL2(), 0.05);
    EXPECT_GT(errors.getH1Seminorm(), 0.2);
  }

  TEST_P(NonlinearPoissonHPTest, QuadratureAndNewtonSensitivity)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const auto a = NonlinearPoissonProblem(mesh).solve<3>();
    const auto b = NonlinearPoissonProblem(mesh, 18).solve<3>(false, 1e-12);
    ASSERT_GT(a.getL2(), 0);
    ASSERT_GT(a.getH1Seminorm(), 0);
    EXPECT_LT(std::abs(b.getL2() / a.getL2() - 1), 1e-6);
    EXPECT_LT(std::abs(b.getH1Seminorm() / a.getH1Seminorm() - 1), 1e-6);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, NonlinearPoissonHPTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
