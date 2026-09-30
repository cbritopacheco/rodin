/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Complex P1 and P2 h-convergence validation for Helmholtz.
 */

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

namespace Rodin::Tests::Convergence::H::Helmholtz
{
  template <class FES, class Exact, class Forcing, class ExactGradient>
  ErrorNorms solveProblem(const LocalMesh& mesh, FES& vh, const Exact& exact,
    const Forcing& forcing, const ExactGradient& exactGradient)
  {
    TrialFunction u(vh);
    TestFunction v(vh);

    auto stiffness = Integral(Grad(u), Grad(v));
    stiffness.setOrder(10);
    auto mass = Integral(u, v);
    mass.setOrder(10);
    auto load = Integral(forcing, v);
    load.setOrder(10);

    Problem helmholtz(u, v);
    helmholtz = stiffness - 0.25 * mass - load + DirichletBC(u, exact);

    CG solver(helmholtz);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());

    return ErrorNorm::compute(mesh, u.getSolution(), exact, exactGradient, 12);
  }

  template <size_t K, class Exact, class Forcing, class ExactGradient>
  ErrorNorms solve(const UniformGridHierarchy& hierarchy, size_t pointsPerAxis,
    const Exact& exact, const Forcing& forcing, const ExactGradient& exactGradient)
  {
    auto mesh = hierarchy.makeMesh(pointsPerAxis);
    if constexpr (K == 1)
    {
      P1<Complex> vh(mesh);
      return solveProblem(mesh, vh, exact, forcing, exactGradient);
    }
    else
    {
      H1<K, Complex> vh(std::integral_constant<size_t, K>{}, mesh);
      return solveProblem(mesh, vh, exact, forcing, exactGradient);
    }
  }

  void expectRates(const ErrorHistory& history, Real minimumL2Rate, Real maximumL2Rate,
    Real minimumH1Rate, Real maximumH1Rate)
  {
    ASSERT_GE(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());

      const auto rates = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "L2 errors " << coarse.getL2() << " -> " << fine.getL2() << ", rate "
        << rates.getL2() << "; H1-seminorm errors " << coarse.getH1Seminorm() << " -> "
        << fine.getH1Seminorm() << ", rate " << rates.getH1Seminorm());
      EXPECT_GT(rates.getL2(), minimumL2Rate);
      EXPECT_LT(rates.getL2(), maximumL2Rate);
      EXPECT_GT(rates.getH1Seminorm(), minimumH1Rate);
      EXPECT_LT(rates.getH1Seminorm(), maximumH1Rate);
    }
  }

  class ComplexHelmholtzHConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  template <size_t K>
  void testPlaneWave(Polytope::Type geometry, std::initializer_list<size_t> levels,
    Real minimumL2Rate, Real maximumL2Rate, Real minimumH1Rate, Real maximumH1Rate)
  {
    const UniformGridHierarchy hierarchy(geometry, levels);
    const size_t dim = hierarchy.getDimension();
    const ComplexFunction exact([dim](const Point& p) {
      Real phase = 0;
      for (size_t i = 0; i < dim; ++i)
        phase += p(i);
      return std::exp(Complex(0, phase));
    });
    const ComplexFunction forcing(
      [dim, &exact](const Point& p) { return (Real(dim) - 0.25) * exact(p); });
    const auto exactGradient = [dim, &exact](const Point& p) {
      Math::SpatialVector<Complex> value(static_cast<std::uint8_t>(dim));
      const Complex derivative = Complex(0, 1) * exact(p);
      for (size_t i = 0; i < dim; ++i)
        value(i) = derivative;
      return value;
    };

    ErrorHistory history;
    for (const size_t level : hierarchy.getLevels())
      history.append(hierarchy.getMeshSize(level),
        solve<K>(hierarchy, level, exact, forcing, exactGradient));
    expectRates(history, minimumL2Rate, maximumL2Rate, minimumH1Rate, maximumH1Rate);
  }

  /** @brief Verifies optimal complex Helmholtz P1 rates. */
  TEST_P(ComplexHelmholtzHConvergenceTest, PlaneWaveHasOptimalP1Rates)
  {
    testPlaneWave<1>(GetParam(), {5, 9, 17}, 1.65, 2.35, 0.75, 1.25);
  }

  /** @brief Verifies optimal complex Helmholtz P2 rates. */
  TEST_P(ComplexHelmholtzHConvergenceTest, PlaneWaveHasOptimalP2Rates)
  {
    testPlaneWave<2>(GetParam(), {3, 5, 9}, 2.45, 3.55, 1.55, 2.45);
  }

  std::string geometryName(const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries, ComplexHelmholtzHConvergenceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    geometryName);
}
