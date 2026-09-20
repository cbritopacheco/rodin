/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/UMFPack.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  using UMFPackLinearSystem =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  TEST(Rodin_Solver_UMFPack, SolvesRawLinearSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    UMFPackLinearSystem system;
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

    Rodin::Solver::UMFPack solver(problem);
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getInfo().factorization, Rodin::Solver::Factorization::Numeric);
    EXPECT_EQ(solver.getInfo().status, 0);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  /// @brief A failed factorization is reported instead of being solved with.
  ///
  /// The solve must not run on a factorization that was never built.
  TEST(Rodin_Solver_UMFPack, ReportsFailedFactorization)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    // Column one is empty, so the matrix is structurally singular. UMFPACK
    // leaves the solution vector untouched in that case, which an unchecked
    // solve would return as if it were the solution.
    UMFPackLinearSystem system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(2, 2) = 4.0;
    system.getVector().resize(3);
    system.getVector() << 1.0, 1.0, 1.0;

    // A sentinel the solver must not overwrite with a bogus solution.
    system.getSolution().resize(3);
    system.getSolution() << 7.0, 7.0, 7.0;

    Rodin::Solver::UMFPack solver(problem);
    solver.solve(system);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(system.getSolution()(0), 7.0);
    EXPECT_EQ(system.getSolution()(1), 7.0);
    EXPECT_EQ(system.getSolution()(2), 7.0);
    // UMFPACK's own status distinguishes the failures, such as a singular
    // matrix from an allocation that could not be made.
    EXPECT_NE(solver.getInfo().status, 0);
    EXPECT_FALSE(solver.getInfo().factorization.has_value());
  }
}
