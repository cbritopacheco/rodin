/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Complex H1 polynomial-degree convergence validation for Helmholtz.
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

namespace Rodin::Tests::Convergence::P::Helmholtz
{
  template <class FES, class Exact, class Forcing, class ExactGradient>
  ErrorNorms solveProblem(const LocalMesh& mesh, FES& space, const Exact& exact,
    const Forcing& forcing, const ExactGradient& exactGradient)
  {
    TrialFunction u(space);
    TestFunction v(space);
    auto stiffness = Integral(Grad(u), Grad(v));
    auto mass = Integral(u, v);
    auto load = Integral(forcing, v);
    stiffness.setOrder(12);
    mass.setOrder(12);
    load.setOrder(12);

    Problem problem(u, v);
    problem = stiffness - 0.25 * mass - load + DirichletBC(u, exact);
    CG solver(problem);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());
    return ErrorNorm::compute(mesh, u.getSolution(), exact, exactGradient, 12);
  }

  template <size_t K, class Exact, class Forcing, class ExactGradient>
  ErrorNorms solve(const LocalMesh& mesh, const Exact& exact, const Forcing& forcing,
    const ExactGradient& exactGradient)
  {
    if constexpr (K == 1)
    {
      P1<Complex> space(mesh);
      return solveProblem(mesh, space, exact, forcing, exactGradient);
    }
    else
    {
      H1<K, Complex> space(std::integral_constant<size_t, K>{}, mesh);
      return solveProblem(mesh, space, exact, forcing, exactGradient);
    }
  }

  class ComplexHelmholtzPConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(ComplexHelmholtzPConvergenceTest, PlaneWaveErrorDecaysWithDegree)
  {
    const UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(2);
    const size_t dim = grid.getDimension();
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
    history.append(1, solve<1>(mesh, exact, forcing, exactGradient))
      .append(2, solve<2>(mesh, exact, forcing, exactGradient))
      .append(3, solve<3>(mesh, exact, forcing, exactGradient))
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
        << history.getSample(i).parameter << "; L2 " << coarse.getL2() << " -> "
        << fine.getL2() << " (rate " << rates.getL2() << "), H1-seminorm "
        << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm() << " (rate "
        << rates.getH1Seminorm() << ')');
      EXPECT_GT(rates.getL2(), 0.25);
      EXPECT_GT(rates.getH1Seminorm(), 0.25);
    }
  }

  std::string geometryName(const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries, ComplexHelmholtzPConvergenceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    geometryName);
}
