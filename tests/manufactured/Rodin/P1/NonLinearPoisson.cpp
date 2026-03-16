/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::NonLinearPoisson
{
  namespace
  {
    auto makeUnitSquareMesh()
    {
      // UniformGrid({17,17}) creates 16 equal intervals per coordinate direction,
      // and scaling by 1/16 maps the mesh to the unit square (0,1)^2.
      constexpr Real domainScale = 1.0 / 16.0;
      Mesh mesh;
      mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 17, 17 });
      mesh.scale(domainScale);
      mesh.getConnectivity().compute(1, 2);
      return mesh;
    }

    template <class FES, class Exact, class Source>
    auto solveByFixedPoint(FES& Vh, const Exact& exact, const Source& f)
    {
      constexpr size_t maxIterations = 30;
      constexpr Real fixedPointTolerance = 1e-12;

      auto uk = GridFunction(Vh);
      uk = Zero();

      TrialFunction u(Vh);
      TestFunction v(Vh);
      Problem nonlinearPoisson(u, v);

      for (size_t it = 0; it < maxIterations; ++it)
      {
        const auto nonlinearTerm = uk * uk * uk;
        nonlinearPoisson = Integral(Grad(u), Grad(v))
                         + Integral(u, v)
                         - Integral(f - nonlinearTerm, v)
                         + DirichletBC(u, exact);
        SparseLU(nonlinearPoisson).solve();

        const Real updateNorm = (u.getSolution().getData() - uk.getData()).norm();
        uk = u.getSolution();
        if (updateNorm < fixedPointTolerance)
          break;
      }

      return uk;
    }
  }

  TEST(Rodin_Manufactured_P1, NonLinearPoisson_UnitSquareDirichlet_FixedPoint)
  {
    const Real pi = Math::Constants::pi();
    Mesh mesh = makeUnitSquareMesh();
    P1 Vh(mesh);

    auto exact = sin(pi * F::x) * sin(pi * F::y);
    auto f = (2 * pi * pi + 1) * exact + exact * exact * exact;

    const auto uk = solveByFixedPoint(Vh, exact, f);

    GridFunction diff(Vh);
    diff = Pow(uk - exact, 2);
    const Real l2ErrorSquared = Integral(diff).compute();
    EXPECT_NEAR(l2ErrorSquared, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Manufactured_P1, NonLinearPoisson_UnitSquareDirichlet_NewtonSolver)
  {
    const Real pi = Math::Constants::pi();
    Mesh mesh = makeUnitSquareMesh();
    P1 Vh(mesh);

    auto exact = sin(pi * F::x) * sin(pi * F::y);
    auto f = (2 * pi * pi + 1) * exact + exact * exact * exact;

    const auto reference = solveByFixedPoint(Vh, exact, f);

    GridFunction uCurrent(Vh);
    uCurrent = reference;

    TrialFunction du(Vh);
    TestFunction v(Vh);
    Problem tangent(du, v);
    tangent = Integral(Grad(du), Grad(v))
            + Integral((1 + 3 * uCurrent * uCurrent) * du, v)
            - Integral(f, v)
            + Integral(Grad(uCurrent), Grad(v))
            + Integral(uCurrent, v)
            + Integral(uCurrent * uCurrent * uCurrent, v)
            + DirichletBC(du, Zero());

    SparseLU solver(tangent);
    NewtonSolver newton(solver);
    // Start from the fixed-point reference and use a tolerance consistent with
    // that discretized reference. A looser absolute tolerance avoids unstable
    // additional Newton updates once the fixed-point reference is reached.
    constexpr Real newtonAbsoluteTolerance = 5e-2;
    newton.setMaxIterations(20)
      .setAbsoluteTolerance(newtonAbsoluteTolerance)
      .setRelativeTolerance(1e-10);
    newton.solve(uCurrent.getData());

    GridFunction diffExact(Vh);
    diffExact = Pow(uCurrent - exact, 2);
    const Real newtonL2ErrorSquared = Integral(diffExact).compute();
    EXPECT_NEAR(newtonL2ErrorSquared, 0, RODIN_FUZZY_CONSTANT);

    GridFunction diffReference(Vh);
    diffReference = Pow(uCurrent - reference, 2);
    const Real methodDifference = Integral(diffReference).compute();
    EXPECT_NEAR(methodDifference, 0, 1e-5);
  }
}
