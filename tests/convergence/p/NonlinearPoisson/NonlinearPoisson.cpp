/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Semilinear Poisson degree-refinement verification. */

#include "../../NonlinearPoisson.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::P::NonlinearPoisson
{
  class NonlinearPoissonPTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(NonlinearPoissonPTest, ErrorsDecayAtEveryDegree)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const NonlinearPoissonProblem problem(mesh);
    const auto e1 = problem.solve<1>(), e2 = problem.solve<2>();
    const auto e3 = problem.solve<3>(), e4 = problem.solve<4>();
    ErrorHistory history;
    history.append(1, e1);
    history.append(2, e2);
    history.append(3, e3);
    history.append(4, e4);
    for (size_t i = 1; i < 4; ++i)
    {
      const auto& a = history.getSample(i - 1).error;
      const auto& b = history.getSample(i).error;
      ASSERT_TRUE(a.isFinite());
      ASSERT_TRUE(b.isFinite());
      ASSERT_GT(b.getL2(), 0);
      ASSERT_GT(b.getH1Seminorm(), 0);
      ASSERT_GT(a.getL2(), b.getL2());
      ASSERT_GT(a.getH1Seminorm(), b.getH1Seminorm());
      SCOPED_TRACE(::testing::Message()
        << "degree=" << i + 1 << " L2 " << a.getL2() << " -> " << b.getL2() << " H1 "
        << a.getH1Seminorm() << " -> " << b.getH1Seminorm());
      EXPECT_GT(std::log(a.getL2() / b.getL2()), 0.1);
      EXPECT_GT(std::log(a.getH1Seminorm() / b.getH1Seminorm()), 0.1);
    }
  }

  TEST_P(NonlinearPoissonPTest, TangentsAgreeWithResidualAndRejectWrongDerivative)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const NonlinearPoissonProblem problem(mesh);
    EXPECT_LT(problem.tangentDefect<1>(), 1e-6);
    EXPECT_LT(problem.tangentDefect<2>(), 1e-6);
    EXPECT_LT(problem.tangentDefect<3>(), 1e-6);
    EXPECT_LT(problem.tangentDefect<4>(), 1e-6);
    EXPECT_GT(problem.tangentDefect<4>(true), 1e-3);
  }

  TEST_P(NonlinearPoissonPTest, RejectsOmittedCubicReaction)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const NonlinearPoissonProblem problem(mesh, 16, 4);
    const auto baseline = problem.solve<4>();
    EXPECT_LT(baseline.getL2(), 0.05);
    EXPECT_LT(baseline.getH1Seminorm(), 0.2);
    const auto errors = problem.solve<4>(true);
    EXPECT_GT(errors.getL2(), 0.05);
    EXPECT_GT(errors.getH1Seminorm(), 0.2);
  }

  TEST_P(NonlinearPoissonPTest, QuadratureAndNewtonSensitivity)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const auto a = NonlinearPoissonProblem(mesh).solve<4>();
    const auto b = NonlinearPoissonProblem(mesh, 18).solve<4>(false, 1e-12);
    ASSERT_GT(a.getL2(), 0);
    ASSERT_GT(a.getH1Seminorm(), 0);
    EXPECT_LT(std::abs(b.getL2() / a.getL2() - 1), 1e-6);
    EXPECT_LT(std::abs(b.getH1Seminorm() / a.getH1Seminorm() - 1), 1e-6);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, NonlinearPoissonPTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
