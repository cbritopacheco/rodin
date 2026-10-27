/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Combined h/p convergence for variable conductivity. */

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

namespace Rodin::Tests::Convergence::HP::Conductivity
{
  template <size_t K, class Exact, class Source, class Gradient>
  ErrorNorms solve(const LocalMesh& mesh, const Exact& exact, const Source& source,
    const Gradient& gradient)
  {
    H1 space(std::integral_constant<size_t, K>{}, mesh);
    TrialFunction u(space);
    TestFunction v(space);
    const size_t dim = mesh.getDimension();
    const RealFunction gamma([dim](const Point& p) {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value += p(i);
      return value;
    });
    auto stiffness = Integral(gamma * Grad(u), Grad(v));
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

  class ConductivityHPConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(ConductivityHPConvergenceTest, AnalyticSolutionImproves)
  {
    const UniformGrid grid(GetParam());
    const size_t dim = grid.getDimension();
    const RealFunction exact([dim](const Point& p) {
      Real exponent = 0;
      for (size_t i = 0; i < dim; ++i)
        exponent += p(i);
      return std::exp(exponent);
    });
    const RealFunction source([dim, &exact](const Point& p) {
      Real gamma = 1;
      for (size_t i = 0; i < dim; ++i)
        gamma += p(i);
      return -Real(dim) * (1 + gamma) * exact(p);
    });
    const VectorFunction gradient(dim, [dim, &exact](const Point& p) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
        value(i) = exact(p);
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
      const auto rates = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "interval " << i << ": L2 " << coarse.getL2() << " -> " << fine.getL2()
        << ", H1 " << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
        << ", ratios " << rates.getL2() << ", " << rates.getH1Seminorm());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      EXPECT_GT(rates.getL2(), 1.9);
      EXPECT_GT(rates.getH1Seminorm(), 0.9);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries, ConductivityHPConvergenceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
