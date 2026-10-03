/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Coupled reaction-diffusion combined mesh/degree validation. */

#include <array>
#include <gtest/gtest.h>

#include "../../ReactionDiffusion.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::HP::ReactionDiffusion
{
  template <size_t K>
  std::array<ErrorNorms, 2> solve(const LocalMesh& mesh,
    ReactionDiffusionData::Field field, bool omitCoupling = false,
    size_t quadratureOrder = 16, Real tolerance = 1e-13)
  {
    const ReactionDiffusionData data(mesh.getDimension(), field);
    const auto exactU = data.getSolution(0), exactW = data.getSolution(1);
    const auto sourceU = data.getSource(0), sourceW = data.getSource(1);
    const auto gradientU = data.getGradient(0), gradientW = data.getGradient(1);
    H1<K, Real> space(std::integral_constant<size_t, K>{}, mesh);
    TrialFunction u(space), w(space);
    TestFunction v(space), z(space);
    const Real alpha = omitCoupling ? 0 : 0.2;
    auto aUU = Integral(Grad(u), Grad(v));
    auto aWW = Integral(2 * Grad(w), Grad(z));
    auto rUU = Integral(u, v);
    auto rWW = Integral(w, z);
    auto rUW = Integral(alpha * w, v);
    auto rWU = Integral(alpha * u, z);
    auto bU = Integral(sourceU, v);
    auto bW = Integral(sourceW, z);
    aUU.setOrder(quadratureOrder);
    aWW.setOrder(quadratureOrder);
    rUU.setOrder(quadratureOrder);
    rWW.setOrder(quadratureOrder);
    rUW.setOrder(quadratureOrder);
    rWU.setOrder(quadratureOrder);
    bU.setOrder(quadratureOrder);
    bW.setOrder(quadratureOrder);
    Problem problem(u, w, v, z);
    problem = aUU + rUU + rUW - bU + aWW + rWW + rWU - bW + DirichletBC(u, exactU) +
      DirichletBC(w, exactW);
    Solver::CG solver(problem);
    solver.setTolerance(tolerance).setMaxIterations(50000).solve();
    EXPECT_TRUE(solver.success());
    return {ErrorNorm::compute(mesh, u.getSolution(), exactU, gradientU, quadratureOrder),
      ErrorNorm::compute(mesh, w.getSolution(), exactW, gradientW, quadratureOrder)};
  }

  class ReactionDiffusionHPTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(ReactionDiffusionHPTest, BothFieldsImproveAlongCombinedPath)
  {
    const UniformGrid grid(GetParam());
    std::array<ErrorHistory, 2> histories;
    const auto append = [&](Real h, const auto& errors) {
      for (size_t component = 0; component < 2; ++component)
        histories[component].append(h, errors[component]);
    };
    append(1, solve<1>(grid.makeMesh(2), ReactionDiffusionData::Field::Smooth));
    append(0.5, solve<2>(grid.makeMesh(3), ReactionDiffusionData::Field::Smooth));
    append(0.25, solve<3>(grid.makeMesh(5), ReactionDiffusionData::Field::Smooth));
    for (size_t component = 0; component < 2; ++component)
    {
      SCOPED_TRACE(::testing::Message() << "component=" << component);
      const auto& history = histories[component];
      ASSERT_EQ(history.getSize(), 3u);
      for (size_t i = 1; i < history.getSize(); ++i)
      {
        const auto& coarse = history.getSample(i - 1).error;
        const auto& fine = history.getSample(i).error;
        ASSERT_TRUE(coarse.isFinite());
        ASSERT_TRUE(fine.isFinite());
        ASSERT_GT(fine.getL2(), 0);
        ASSERT_GT(fine.getH1Seminorm(), 0);
        ASSERT_GT(coarse.getL2(), fine.getL2());
        ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
        const auto rate = history.getAlgebraicRates(i);
        SCOPED_TRACE(::testing::Message()
          << "interval=" << i << " L2 " << coarse.getL2() << " -> " << fine.getL2()
          << " H1 " << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
          << " effective rates " << rate.getL2() << ", " << rate.getH1Seminorm());
        EXPECT_GT(rate.getL2(), 1.9);
        EXPECT_GT(rate.getH1Seminorm(), 0.9);
      }
    }
  }

  TEST_P(ReactionDiffusionHPTest, P3ReproducesQuadraticPatch)
  {
    const UniformGrid grid(GetParam());
    for (const auto& error :
      solve<3>(grid.makeMesh(3), ReactionDiffusionData::Field::Quadratic))
    {
      EXPECT_LT(error.getL2(), 1e-9);
      EXPECT_LT(error.getH1Seminorm(), 1e-9);
    }
  }

  TEST_P(ReactionDiffusionHPTest, P3PatchRejectsOmittedCoupling)
  {
    const UniformGrid grid(GetParam());
    for (const auto& error :
      solve<3>(grid.makeMesh(3), ReactionDiffusionData::Field::Affine, true))
    {
      EXPECT_GT(error.getL2(), 1e-3);
      EXPECT_GT(error.getH1Seminorm(), 1e-2);
    }
  }

  TEST_P(ReactionDiffusionHPTest, P3QuadratureAndSolverSensitivity)
  {
    const UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(3);
    const auto baseline = solve<3>(mesh, ReactionDiffusionData::Field::Smooth);
    const auto refined =
      solve<3>(mesh, ReactionDiffusionData::Field::Smooth, false, 18, 1e-14);
    for (size_t component = 0; component < 2; ++component)
    {
      ASSERT_GT(baseline[component].getL2(), 0);
      ASSERT_GT(baseline[component].getH1Seminorm(), 0);
      EXPECT_LT(
        std::abs(refined[component].getL2() / baseline[component].getL2() - 1), 1e-6);
      EXPECT_LT(
        std::abs(
          refined[component].getH1Seminorm() / baseline[component].getH1Seminorm() - 1),
        1e-6);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, ReactionDiffusionHPTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
