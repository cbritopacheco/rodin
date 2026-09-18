/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/ParU.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  TEST(Rodin_Solver_ParU, SolvesRawLinearSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    Rodin::Solver::ParU solver(problem);
    using LinearSystemType =
      Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    LinearSystemType system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 4.0;
    system.getOperator().insert(1, 0) = 2.0;
    system.getOperator().insert(0, 1) = 1.0;
    system.getOperator().insert(1, 1) = 5.0;
    system.getOperator().insert(2, 1) = 3.0;
    system.getOperator().insert(1, 2) = 1.0;
    system.getOperator().insert(2, 2) = 6.0;
    system.getVector().resize(3);
    system.getVector() << 6.0, 15.0, 24.0;

    solver
      .setMaxThreads(1)
      .setOrdering(decltype(solver)::Ordering::AMD)
      .solve(system);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  TEST(Rodin_Solver_ParU, SolvesVariationalProblemWithCTAD)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    RealFunction f = 1.0;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) +
      DirichletBC(u, Zero());

    Rodin::Solver::ParU solver(poisson);
    static_assert(std::is_same_v<
      typename FormLanguage::Traits<decltype(solver)>::LinearSystemType,
      typename decltype(poisson)::LinearSystemType>);
    solver
      .setMaxThreads(1)
      .setOrdering(decltype(solver)::Ordering::AMD)
      .solve();

    const auto& system = poisson.getLinearSystem();
    EXPECT_LT(
      (system.getOperator() * system.getSolution() - system.getVector()).norm(),
      1e-11);
  }
}
