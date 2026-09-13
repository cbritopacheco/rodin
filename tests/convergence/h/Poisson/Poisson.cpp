/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief P1 h-convergence validation for the Poisson equation.
 *
 * On a shape-regular mesh family and for a smooth exact solution, the
 * conforming P1 Galerkin approximation satisfies
 * @f[
 *   \lVert u-u_h\rVert_{H^1(\Omega)}=O(h), \qquad
 *   \lVert u-u_h\rVert_{L^2(\Omega)}=O(h^2).
 * @f]
 * The tests below measure both errors by independent order-eight quadrature,
 * require strict error reduction, and bound every observed refinement rate.
 * They cover every positive-dimensional cell geometry supported by
 * Geometry::Mesh::UniformGrid.
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

namespace Rodin::Tests::Convergence::H::Poisson
{
  template <class Exact, class Forcing, class ExactGradient>
  ErrorNorms solve(
    const UniformGridHierarchy& hierarchy,
    size_t pointsPerAxis,
    const Exact& exact,
    const Forcing& forcing,
    const ExactGradient& exactGradient)
  {
    auto mesh = hierarchy.makeMesh(pointsPerAxis);
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    auto load = Integral(forcing, v);
    load.setOrder(8);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - load
            + DirichletBC(u, exact);

    CG solver(poisson);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());

    return ErrorNorm::compute(mesh, u.getSolution(), exact, exactGradient);
  }

  void expectRate(
    const ErrorHistory& history,
    Real minimumL2Rate,
    Real maximumL2Rate,
    Real minimumH1Rate,
    Real maximumH1Rate)
  {
    ASSERT_GE(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      SCOPED_TRACE(::testing::Message() << "refinement interval " << i);
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(fine.getL2(), 0);
      ASSERT_GT(fine.getH1Seminorm(), 0);
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());

      const auto rates = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "L2 errors " << coarse.getL2() << " -> " << fine.getL2()
        << ", rate " << rates.getL2() << "; H1 seminorm errors "
        << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
        << ", rate " << rates.getH1Seminorm());
      EXPECT_GT(rates.getL2(), minimumL2Rate);
      EXPECT_LT(rates.getL2(), maximumL2Rate);
      EXPECT_GT(rates.getH1Seminorm(), minimumH1Rate);
      EXPECT_LT(rates.getH1Seminorm(), maximumH1Rate);
    }
  }

  class PoissonHConvergenceTest
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * @brief Verifies the P1 affine patch test on every UniformGrid geometry.
   *
   * The exact solution @f$u=1+\sum_i x_i@f$ belongs to the discrete space and
   * satisfies @f$-\Delta u=0@f$. Galerkin consistency therefore requires
   * roundoff-level value and gradient errors, independently of @f$h@f$.
   */
  TEST_P(PoissonHConvergenceTest, AffineSolutionIsExact)
  {
    const auto geometry = GetParam();
    const UniformGridHierarchy hierarchy(geometry, {5, 9, 17});
    const size_t dim = hierarchy.getDimension();
    const RealFunction exact([dim](const Point& p)
      {
        Real value = 1;
        for (size_t i = 0; i < dim; ++i)
          value += p(i);
        return value;
      });
    const RealFunction forcing(0.0);
    const VectorFunction exactGradient(dim, [dim](const Point&)
      {
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        for (size_t i = 0; i < dim; ++i)
          value(i) = 1;
        return value;
      });

    const auto error = solve(hierarchy, dim < 3 ? 9 : 5,
      exact, forcing, exactGradient);
    EXPECT_LT(error.getL2(), 1e-10);
    EXPECT_LT(error.getH1Seminorm(), 1e-10);
  }

  /**
   * @brief Verifies optimal P1 rates for a homogeneous Dirichlet solution.
   *
   * For @f$u=\prod_i\sin(\pi x_i)@f$, the load is
   * @f$f=d\pi^2u@f$. The solution is smooth and vanishes on the whole boundary,
   * so the expected rates are two in L2 and one in the H1 seminorm.
   */
  TEST_P(PoissonHConvergenceTest, HomogeneousDirichletHasOptimalRates)
  {
    const auto geometry = GetParam();
    const UniformGridHierarchy hierarchy(geometry, {5, 9, 17});
    const size_t dim = hierarchy.getDimension();
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
    for (const size_t level : hierarchy.getLevels())
      history.append(hierarchy.getMeshSize(level),
        solve(hierarchy, level, exact, forcing, exactGradient));

    expectRate(history, 1.65, 2.35, 0.75, 1.25);
  }

  /**
   * @brief Verifies optimal P1 rates with nonhomogeneous Dirichlet data.
   *
   * For @f$u=\exp(\sum_i x_i)@f$, @f$f=-d u@f$ and the exact trace is
   * imposed on the complete boundary. This separates boundary elimination
   * errors from the zero-trace case while preserving full regularity.
   */
  TEST_P(PoissonHConvergenceTest, NonhomogeneousDirichletHasOptimalRates)
  {
    const auto geometry = GetParam();
    const UniformGridHierarchy hierarchy(geometry, {5, 9, 17});
    const size_t dim = hierarchy.getDimension();
    const RealFunction exact([dim](const Point& p)
      {
        Real sum = 0;
        for (size_t i = 0; i < dim; ++i)
          sum += p(i);
        return std::exp(sum);
      });
    const RealFunction forcing([dim](const Point& p)
      {
        Real sum = 0;
        for (size_t i = 0; i < dim; ++i)
          sum += p(i);
        return -Real(dim) * std::exp(sum);
      });
    const VectorFunction exactGradient(dim, [dim](const Point& p)
      {
        Real sum = 0;
        for (size_t i = 0; i < dim; ++i)
          sum += p(i);
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        for (size_t i = 0; i < dim; ++i)
          value(i) = std::exp(sum);
        return value;
      });

    ErrorHistory history;
    for (const size_t level : hierarchy.getLevels())
      history.append(hierarchy.getMeshSize(level),
        solve(hierarchy, level, exact, forcing, exactGradient));

    expectRate(history, 1.65, 2.35, 0.75, 1.25);
  }

  std::string geometryName(
    const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGridHierarchy::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(
    AllUniformGridGeometries,
    PoissonHConvergenceTest,
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
