/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <limits>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::NonLinearPoisson
{
  TEST(Rodin_Manufactured_P1, NonLinearPoisson_UnitSquareDirichlet)
  {
    const Real pi = Math::Constants::pi();
    // Fixed-point iteration is contractive for this smooth manufactured setup
    // and converges well within this cap.
    constexpr size_t maxIterations = 30;
    // Tight tolerance keeps the nonlinear iteration error below FE discretization error.
    constexpr Real fixedPointTolerance = 1e-12;
    // UniformGrid({17,17}) generates 17 nodes per axis; scale by 1/16 to map to (0,1)^2.
    constexpr Real domainScale = 1.0 / 16.0;

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 17, 17 });
    mesh.scale(domainScale);
    mesh.getConnectivity().compute(1, 2);

    P1 Vh(mesh);

    auto exact = sin(pi * F::x) * sin(pi * F::y);
    auto f = (2 * pi * pi + 1) * exact + exact * exact * exact;

    GridFunction uk(Vh);
    uk = Zero();

    TrialFunction u(Vh);
    TestFunction v(Vh);
    Problem nonlinearPoisson(u, v);

    bool converged = false;
    size_t iterationCount = 0;
    Real lastUpdateNorm = std::numeric_limits<Real>::infinity();
    for (size_t it = 0; it < maxIterations; ++it)
    {
      iterationCount = it + 1;
      const auto nonlinearTerm = uk * uk * uk;
      nonlinearPoisson = Integral(Grad(u), Grad(v))
                       + Integral(u, v)
                       - Integral(f - nonlinearTerm, v)
                       + DirichletBC(u, exact);
      SparseLU(nonlinearPoisson).solve();

      const Real updateNorm = (u.getSolution().getData() - uk.getData()).norm();
      lastUpdateNorm = updateNorm;
      uk = u.getSolution();
      if (updateNorm < fixedPointTolerance)
      {
        converged = true;
        break;
      }
    }
    EXPECT_TRUE(converged)
      << "Fixed-point solver failed to converge in " << iterationCount
      << " iterations (last update norm = " << lastUpdateNorm << ").";

    GridFunction diff(Vh);
    diff = Pow(uk - exact, 2);
    const Real l2ErrorSquared = Integral(diff).compute();
    EXPECT_NEAR(l2ErrorSquared, 0, RODIN_FUZZY_CONSTANT);
  }
}
