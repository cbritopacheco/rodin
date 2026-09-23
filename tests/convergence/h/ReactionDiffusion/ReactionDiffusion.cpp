/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Coupled reaction–diffusion rates on all cell geometries. */

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::ReactionDiffusion
{
  constexpr Real alpha = 0.2;

  template <size_t K, class ExactU, class ExactW, class SourceU,
    class SourceW, class GradientU, class GradientW>
  std::pair<ErrorNorms, ErrorNorms> solve(const LocalMesh& mesh,
    const ExactU& exactU, const ExactW& exactW,
    const SourceU& sourceU, const SourceW& sourceW,
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
    aUU.setOrder(12);
    aWW.setOrder(12);
    rUU.setOrder(12);
    rWW.setOrder(12);
    rUW.setOrder(12);
    rWU.setOrder(12);
    bU.setOrder(12);
    bW.setOrder(12);
    Problem problem(u, w, v, z);
    problem = aUU + rUU + rUW - bU
            + aWW + rWW + rWU - bW
            + DirichletBC(u, exactU)
            + DirichletBC(w, exactW);
    CG solver(problem);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());
    return {
      ErrorNorm::compute(mesh, u.getSolution(), exactU, gradientU, 12),
      ErrorNorm::compute(mesh, w.getSolution(), exactW, gradientW, 12)
    };
  }

  template <size_t K>
  void checkRates(Polytope::Type geometry)
  {
    UniformGridHierarchy hierarchy(geometry,
      K == 1 ? std::initializer_list<size_t>{5, 9, 17}
             : std::initializer_list<size_t>{3, 5, 9});
    const size_t dim = hierarchy.getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exactU([dim, pi](const Point& p)
    {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const RealFunction exactW([dim, pi](const Point& p)
    {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value *= std::cos(pi * p(i));
      return value;
    });
    const VectorFunction gradientU(dim, [dim, pi](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
      {
        value(i) = pi * std::cos(pi * p(i));
        for (size_t j = 0; j < dim; ++j)
          if (j != i)
            value(i) *= std::sin(pi * p(j));
      }
      return value;
    });
    const VectorFunction gradientW(dim, [dim, pi](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
      {
        value(i) = -pi * std::sin(pi * p(i));
        for (size_t j = 0; j < dim; ++j)
          if (j != i)
            value(i) *= std::cos(pi * p(j));
      }
      return value;
    });
    const RealFunction sourceU([&](const Point& p)
    {
      return (Real(dim) * pi * pi + 1) * exactU(p) + alpha * exactW(p);
    });
    const RealFunction sourceW([&](const Point& p)
    {
      return (Real(dim) * pi * pi + 1) * exactW(p) + alpha * exactU(p);
    });
    ErrorHistory uHistory, wHistory;
    for (size_t level : hierarchy.getLevels())
    {
      const auto mesh = hierarchy.makeMesh(level);
      const auto [uError, wError] = solve<K>(mesh,
        exactU, exactW, sourceU, sourceW, gradientU, gradientW);
      const Real h = hierarchy.getMeshSize(level);
      uHistory.append(h, uError);
      wHistory.append(h, wError);
    }
    for (const auto* history : {&uHistory, &wHistory})
    {
      for (size_t i = 1; i < history->getSize(); ++i)
      {
        const auto& coarse = history->getSample(i - 1).error;
        const auto& fine = history->getSample(i).error;
        ASSERT_TRUE(coarse.isFinite());
        ASSERT_TRUE(fine.isFinite());
        ASSERT_GT(coarse.getL2(), fine.getL2());
        ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
        const auto rate = history->getAlgebraicRates(i);
        SCOPED_TRACE(::testing::Message() << "L2 " << coarse.getL2()
          << " -> " << fine.getL2() << ", H1 " << coarse.getH1Seminorm()
          << " -> " << fine.getH1Seminorm() << ", rates " << rate.getL2()
          << ", " << rate.getH1Seminorm());
        EXPECT_GT(rate.getL2(), Real(K) + 0.5);
        EXPECT_LT(rate.getL2(), Real(K) + 1.5);
        EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.3);
        EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.5);
      }
    }
  }

  class ReactionDiffusionTest : public ::testing::TestWithParam<Polytope::Type> {};
  TEST_P(ReactionDiffusionTest, P1OptimalRates) { checkRates<1>(GetParam()); }
  TEST_P(ReactionDiffusionTest, P2OptimalRates) { checkRates<2>(GetParam()); }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, ReactionDiffusionTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
