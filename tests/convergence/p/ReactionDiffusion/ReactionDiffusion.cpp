/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Coupled reaction--diffusion p-rates on every cell geometry. */

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <string>
#include <utility>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::P::ReactionDiffusion
{
  constexpr Real alpha = 0.2;

  template <size_t K, class ExactU, class ExactW, class SourceU, class SourceW,
    class GradientU, class GradientW>
  std::pair<ErrorNorms, ErrorNorms> solve(const LocalMesh& mesh, const ExactU& exactU,
    const ExactW& exactW, const SourceU& sourceU, const SourceW& sourceW,
    const GradientU& gradientU, const GradientW& gradientW)
  {
    H1 space(std::integral_constant<size_t, K>{}, mesh);
    TrialFunction u(space), w(space);
    TestFunction v(space), z(space);
    auto aUU = Integral(Grad(u), Grad(v));
    auto aWW = Integral(Grad(w), Grad(z));
    auto rUU = Integral(u, v);
    auto rWW = Integral(w, z);
    auto rUW = Integral(alpha * w, v);
    auto rWU = Integral(alpha * u, z);
    auto bU = Integral(sourceU, v);
    auto bW = Integral(sourceW, z);
    aUU.setOrder(16);
    aWW.setOrder(16);
    rUU.setOrder(16);
    rWW.setOrder(16);
    rUW.setOrder(16);
    rWU.setOrder(16);
    bU.setOrder(16);
    bW.setOrder(16);
    Problem problem(u, w, v, z);
    problem = aUU + rUU + rUW - bU + aWW + rWW + rWU - bW + DirichletBC(u, exactU) +
      DirichletBC(w, exactW);
    CG solver(problem);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());
    return {ErrorNorm::compute(mesh, u.getSolution(), exactU, gradientU, 16),
      ErrorNorm::compute(mesh, w.getSolution(), exactW, gradientW, 16)};
  }

  class ReactionDiffusionTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(ReactionDiffusionTest, BothFieldsConvergeWithDegree)
  {
    const UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(2);
    const size_t dim = grid.getDimension();
    const RealFunction exactU([dim](const Point& p) {
      Real sum = 0;
      for (size_t i = 0; i < dim; ++i)
        sum += p(i);
      return std::exp(sum);
    });
    const RealFunction exactW([dim](const Point& p) {
      Real sum = 0;
      for (size_t i = 0; i < dim; ++i)
        sum += p(i);
      return std::exp(-sum);
    });
    const VectorFunction gradientU(dim, [dim, &exactU](const Point& p) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
        value(i) = exactU(p);
      return value;
    });
    const VectorFunction gradientW(dim, [dim, &exactW](const Point& p) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
        value(i) = -exactW(p);
      return value;
    });
    const RealFunction sourceU([dim, &exactU, &exactW](const Point& p) {
      return (Real(1) - Real(dim)) * exactU(p) + alpha * exactW(p);
    });
    const RealFunction sourceW([dim, &exactU, &exactW](const Point& p) {
      return (Real(1) - Real(dim)) * exactW(p) + alpha * exactU(p);
    });

    ErrorHistory uHistory, wHistory;
    const auto append = [&](size_t degree, const auto& errors) {
      uHistory.append(Real(degree), errors.first);
      wHistory.append(Real(degree), errors.second);
    };
    append(1, solve<1>(mesh, exactU, exactW, sourceU, sourceW, gradientU, gradientW));
    append(2, solve<2>(mesh, exactU, exactW, sourceU, sourceW, gradientU, gradientW));
    append(3, solve<3>(mesh, exactU, exactW, sourceU, sourceW, gradientU, gradientW));
    append(4, solve<4>(mesh, exactU, exactW, sourceU, sourceW, gradientU, gradientW));

    for (const auto* history : {&uHistory, &wHistory})
    {
      ASSERT_EQ(history->getSize(), 4);
      for (size_t i = 1; i < history->getSize(); ++i)
      {
        const auto& coarse = history->getSample(i - 1).error;
        const auto& fine = history->getSample(i).error;
        ASSERT_TRUE(coarse.isFinite());
        ASSERT_TRUE(fine.isFinite());
        const auto rate = history->getExponentialRates(i);
        SCOPED_TRACE(::testing::Message()
          << "degrees " << i << " -> " << i + 1 << "; L2 " << coarse.getL2() << " -> "
          << fine.getL2() << " (rate " << rate.getL2() << "), H1 "
          << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm() << " (rate "
          << rate.getH1Seminorm() << ')');
        EXPECT_GT(rate.getL2(), 0.25);
        EXPECT_GT(rate.getH1Seminorm(), 0.25);
      }
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, ReactionDiffusionTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
