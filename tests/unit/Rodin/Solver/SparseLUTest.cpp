/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  using SparseLULinearSystem =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  TEST(Rodin_Solver_SparseLU, SolvesRawLinearSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    SparseLULinearSystem system;
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

    Rodin::Solver::SparseLU solver(problem);
    solver.solve(system);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  /// @brief A failed factorization is reported instead of being solved with.
  TEST(Rodin_Solver_SparseLU, ReportsFailedFactorization)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    // Column one is empty, so the matrix is structurally singular.
    SparseLULinearSystem system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(2, 2) = 4.0;
    system.getVector().resize(3);
    system.getVector() << 1.0, 1.0, 1.0;

    Rodin::Solver::SparseLU solver(problem);
    EXPECT_ANY_THROW(solver.solve(system));
  }
}
