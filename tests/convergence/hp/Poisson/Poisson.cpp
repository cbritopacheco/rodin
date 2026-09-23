/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Combined h/p convergence for analytic Poisson data. */

#include <cmath>
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

namespace Rodin::Tests::Convergence::HP::Poisson
{
  template <size_t K, class Exact, class Source, class Gradient>
  ErrorNorms solve(const LocalMesh& mesh, const Exact& exact,
    const Source& source, const Gradient& gradient)
  {
    H1 space(std::integral_constant<size_t, K>{}, mesh);
    TrialFunction u(space);
    TestFunction v(space);
    auto stiffness = Integral(Grad(u), Grad(v));
    auto load = Integral(source, v);
    stiffness.setOrder(12);
    load.setOrder(12);
    Problem problem(u, v);
    problem = stiffness - load + DirichletBC(u, exact);
    CG solver(problem);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());
    return ErrorNorm::compute(mesh, u.getSolution(), exact, gradient, 12);
  }

  class HPPoissonTest : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(HPPoissonTest, AnalyticSolutionImprovesUnderCombinedRefinement)
  {
    UniformGrid grid(GetParam());
    const size_t dim = grid.getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
    {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const RealFunction source([dim, pi](const Point& p)
    {
      Real value = Real(dim) * pi * pi;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const VectorFunction gradient(dim, [dim, pi](const Point& p)
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
    history.append(1, solve<1>(grid.makeMesh(2), exact, source, gradient))
           .append(0.5, solve<2>(grid.makeMesh(3), exact, source, gradient))
           .append(0.25, solve<3>(grid.makeMesh(5), exact, source, gradient));
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message() << "interval " << i
        << ": L2 " << coarse.getL2() << " -> " << fine.getL2()
        << ", H1 " << coarse.getH1Seminorm() << " -> "
        << fine.getH1Seminorm() << ", ratios " << rate.getL2()
        << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), 1.9);
      EXPECT_GT(rate.getH1Seminorm(), 0.9);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, HPPoissonTest,
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
