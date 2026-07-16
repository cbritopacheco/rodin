/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using Rodin::Solver::CG;
using Rodin::Solver::SparseLU;

namespace Rodin::Tests::Unit::InitialGuess
{
  /**
   * A full assemble() must seed the linear system's solution vector (the
   * initial guess) from the trial function data.
   */
  TEST(Rodin_Solver_InitialGuess, AssembleSeedsGuessFromTrialFunction)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    RealFunction f = 1.0;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    u.getSolution() = [](const Geometry::Point& p)
    { return p.x() + 2.0 * p.y(); };

    poisson.assemble();

    const auto& guess = poisson.getLinearSystem().getSolution();
    const auto& data = u.getSolution().getData();
    ASSERT_EQ(guess.size(), data.size());
    EXPECT_EQ((guess - data).norm(), 0.0);
  }

  /**
   * An iterative solver must use the seeded guess: with the exact solution
   * as the guess, CG must succeed within an iteration budget that is
   * insufficient from a zero guess.
   */
  TEST(Rodin_Solver_InitialGuess, IterativeSolverHonorsGuess)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    RealFunction f = 1.0;

    auto define = [&](auto& pb)
    {
      pb = Integral(Grad(u), Grad(v))
         - Integral(f, v)
         + DirichletBC(u, Zero());
    };

    // Reference solution with a direct solver.
    Problem direct(u, v);
    define(direct);
    SparseLU lu(direct);
    direct.solve(lu);
    const Math::Vector<Real> exact = u.getSolution().getData();
    ASSERT_GT(exact.norm(), 0.0);

    // u's solution holds the exact data; assemble seeds it as the guess and
    // a single CG iteration suffices.
    Problem warm(u, v);
    define(warm);
    CG cg(warm);
    cg.setTolerance(1e-10).setMaxIterations(1);
    warm.solve(cg);
    EXPECT_TRUE(cg.success());
    EXPECT_NEAR((u.getSolution().getData() - exact).norm(), 0.0, 1e-10);

    // From a zero guess, the same iteration budget must fail.
    Problem cold(u, v);
    define(cold);
    u.getSolution() = Zero();
    CG coldCG(cold);
    coldCG.setTolerance(1e-10).setMaxIterations(1);
    cold.solve(coldCG);
    EXPECT_FALSE(coldCG.success());
  }

  /**
   * The Newton increment solve always starts from a zero guess, regardless
   * of what the assembly seeded: a direct and an iterative solver must
   * converge to the same Newton solution.
   */
  TEST(Rodin_Solver_InitialGuess, NewtonIncrementStartsFromZero)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 8, 8 });
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);

    auto solveWith = [&](auto makeSolver) -> Math::Vector<Real>
    {
      TrialFunction du(vh);
      TestFunction v(vh);

      GridFunction u(vh);
      u = Zero();

      RealFunction f = 1.0;

      // Newton form of the linear Poisson problem: J du = -F(u).
      Problem newton(du, v);
      newton = Integral(Grad(du), Grad(v))
             + Integral(Grad(u), Grad(v))
             - Integral(f, v)
             + DirichletBC(du, Zero());

      auto linearSolver = makeSolver(newton);
      Solver::NewtonSolver solver(linearSolver);
      solver.setMaxIterations(10)
            .setAbsoluteTolerance(1e-12)
            .setRelativeTolerance(1e-10);
      solver.solve(u);
      EXPECT_TRUE(solver.getReport().converged);
      return u.getData();
    };

    const auto direct =
      solveWith([](auto& pb) { return SparseLU(pb); });

    const auto iterative =
      solveWith([](auto& pb)
      {
        CG cg(pb);
        cg.setTolerance(1e-12).setMaxIterations(2000);
        return cg;
      });

    EXPECT_GT(direct.norm(), 0.0);
    EXPECT_NEAR((direct - iterative).norm(), 0.0, 1e-8);
  }
}
