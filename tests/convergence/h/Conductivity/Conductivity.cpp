/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief All-geometry h-rates for variable-coefficient conductivity. */

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

namespace Rodin::Tests::Convergence::H::Conductivity
{
  template <size_t K, class Exact, class Source, class Gradient>
  ErrorNorms solve(const LocalMesh& mesh, const Exact& exact,
    const Source& source, const Gradient& gradient)
  {
    H1 space(std::integral_constant<size_t, K>{}, mesh);
    TrialFunction u(space);
    TestFunction v(space);
    const size_t dim = mesh.getDimension();
    const RealFunction gamma([dim](const Point& p)
    {
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

  class ConductivityTest : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(ConductivityTest, AffinePatchIsExact)
  {
    UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(5);
    const size_t dim = grid.getDimension();
    const RealFunction exact([dim](const Point& p)
    {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value += p(i);
      return value;
    });
    const RealFunction source(-Real(dim));
    const VectorFunction gradient(dim, [dim](const Point&)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
        value(i) = 1;
      return value;
    });
    const auto error = solve<1>(mesh, exact, source, gradient);
    EXPECT_LT(error.getL2(), 1e-10);
    EXPECT_LT(error.getH1Seminorm(), 1e-10);
  }

  template <size_t K>
  void checkSmoothRates(Polytope::Type geometry)
  {
    UniformGridHierarchy hierarchy(geometry,
      K == 1 ? std::initializer_list<size_t>{5, 9, 17}
             : std::initializer_list<size_t>{3, 5, 9});
    const size_t dim = hierarchy.getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
    {
      Real value = 1;
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
    const RealFunction source([dim, pi](const Point& p)
    {
      Real gamma = 1;
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
      {
        gamma += p(i);
        value *= std::sin(pi * p(i));
      }
      Real result = Real(dim) * pi * pi * gamma * value;
      for (size_t i = 0; i < dim; ++i)
      {
        Real derivative = pi * std::cos(pi * p(i));
        for (size_t j = 0; j < dim; ++j)
          if (j != i)
            derivative *= std::sin(pi * p(j));
        result -= derivative;
      }
      return result;
    });
    ErrorHistory history;
    for (size_t level : hierarchy.getLevels())
    {
      const auto mesh = hierarchy.makeMesh(level);
      history.append(hierarchy.getMeshSize(level),
        solve<K>(mesh, exact, source, gradient));
    }
    ASSERT_EQ(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message() << "L2 " << coarse.getL2()
        << " -> " << fine.getL2() << ", H1 " << coarse.getH1Seminorm()
        << " -> " << fine.getH1Seminorm() << ", rates " << rate.getL2()
        << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), Real(K) + 0.6);
      EXPECT_LT(rate.getL2(), Real(K) + 1.4);
      EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.25);
      EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.4);
    }
  }

  TEST_P(ConductivityTest, P1OptimalRates) { checkSmoothRates<1>(GetParam()); }
  TEST_P(ConductivityTest, P2OptimalRates) { checkSmoothRates<2>(GetParam()); }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, ConductivityTest,
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
