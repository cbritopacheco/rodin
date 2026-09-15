/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief H1 polynomial-degree convergence validation for Poisson.
 */

#include <cstdint>
#include <string>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::P::Poisson
{
  template <size_t K, class Exact, class Forcing, class ExactGradient>
  ErrorNorms solve(
    const LocalMesh& mesh,
    const Exact& exact,
    const Forcing& forcing,
    const ExactGradient& exactGradient)
  {
    H1 vh(std::integral_constant<size_t, K>{}, mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    auto load = Integral(forcing, v);
    load.setOrder(12);
    auto bilinear = Integral(Grad(u), Grad(v));
    bilinear.setOrder(12);
    Problem poisson(u, v);
    poisson = bilinear
            - load
            + DirichletBC(u, exact);

    CG solver(poisson);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());
    return ErrorNorm::compute(
      mesh, u.getSolution(), exact, exactGradient, 12);
  }

  class PoissonPConvergenceTest
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PoissonPConvergenceTest, QuadraticSolutionIsExactFromDegreeTwo)
  {
    const UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(2);
    const size_t dim = grid.getDimension();
    const RealFunction exact([dim](const Point& p)
      {
        Real value = 0;
        for (size_t i = 0; i < dim; ++i)
          value += p(i) * p(i);
        return value;
      });
    const RealFunction forcing(-2 * Real(dim));
    const VectorFunction exactGradient(dim, [dim](const Point& p)
      {
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        for (size_t i = 0; i < dim; ++i)
          value(i) = 2 * p(i);
        return value;
      });

    const auto p1 = solve<1>(mesh, exact, forcing, exactGradient);
    const auto p2 = solve<2>(mesh, exact, forcing, exactGradient);
    EXPECT_GT(p1.getL2(), 1e-3);
    EXPECT_GT(p1.getH1Seminorm(), 1e-2);
    EXPECT_LT(p2.getL2(), 1e-10);
    EXPECT_LT(p2.getH1Seminorm(), 1e-10);
  }

  /**
   * @brief Verifies exponential p-decay between degrees one and four.
   *
   * For the analytic sine-product solution, approximation theory gives
   * @f$e_p\le C\exp(-\alpha p)@f$. The measured L2 and H1-seminorm decay
   * constants must both be positive on every cell geometry.
   */
  TEST_P(PoissonPConvergenceTest, AnalyticSolutionConvergesExponentially)
  {
    const UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(2);
    const size_t dim = grid.getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
      {
        Real value = 1;
        for (size_t i = 0; i < dim; ++i)
          value *= std::sin(pi * p(i));
        return value;
      });
    const RealFunction forcing([dim, pi](const Point& p)
      {
        Real value = Real(dim) * pi * pi;
        for (size_t i = 0; i < dim; ++i)
          value *= std::sin(pi * p(i));
        return value;
      });
    const VectorFunction exactGradient(dim, [dim, pi](const Point& p)
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

    ErrorHistory history;
    history.append(1, solve<1>(mesh, exact, forcing, exactGradient))
           .append(4, solve<4>(mesh, exact, forcing, exactGradient));

    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      const auto rates = history.getExponentialRates(i);
      SCOPED_TRACE(::testing::Message()
        << "degrees " << history.getSample(i - 1).parameter << " -> "
        << history.getSample(i).parameter << "; L2 " << coarse.getL2()
        << " -> " << fine.getL2() << " (rate " << rates.getL2()
        << "), H1-seminorm " << coarse.getH1Seminorm() << " -> "
        << fine.getH1Seminorm() << " (rate " << rates.getH1Seminorm()
        << ')');
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      EXPECT_GT(rates.getL2(), 0.1);
      EXPECT_GT(rates.getH1Seminorm(), 0.1);
    }
  }

  std::string geometryName(
    const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(
    AllUniformGridGeometries,
    PoissonPConvergenceTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    geometryName);
}
