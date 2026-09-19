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
  static Mesh<Context::Local> makeMesh(Polytope::Type geometry)
  {
    switch (geometry)
    {
      case Polytope::Type::Segment:
      {
        auto mesh = Mesh<Context::Local>::UniformGrid(geometry, {4});
        mesh.getConnectivity().compute(0, 1);
        mesh.getConnectivity().compute(1, 0);
        return mesh;
      }
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Wedge:
      {
        auto mesh = Mesh<Context::Local>::UniformGrid(geometry, {3, 3, 3});
        mesh.getConnectivity().compute(2, 3);
        mesh.getConnectivity().compute(3, 0);
        return mesh;
      }
      default:
      {
        auto mesh = Mesh<Context::Local>::UniformGrid(geometry, {4, 4});
        mesh.getConnectivity().compute(1, 2);
        mesh.getConnectivity().compute(2, 0);
        return mesh;
      }
    }
  }

  TEST(Rodin_Solver_ParU, SharesValuesAndUses32BitIndices)
  {
    Math::SparseMatrix<Real> matrix(2, 2);
    matrix.insert(0, 0) = 2.0;
    matrix.insert(1, 1) = 3.0;
    matrix.makeCompressed();

    cholmod_sparse view = Eigen::viewAsCholmod(matrix);

    EXPECT_EQ(view.itype, CHOLMOD_INT);
    EXPECT_EQ(view.p, static_cast<void*>(matrix.outerIndexPtr()));
    EXPECT_EQ(view.i, static_cast<void*>(matrix.innerIndexPtr()));
    EXPECT_EQ(view.x, static_cast<void*>(matrix.valuePtr()));
  }

  TEST(Rodin_Solver_ParU, SolvesRawLinearSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    Rodin::Solver::ParU solver(problem);
    using enum decltype(solver)::Factorization;
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

    solver.setMaxThreads(1)
      .setOrdering(decltype(solver)::Ordering::AMD)
      .factorize(system);
    EXPECT_TRUE(solver.hasFactorization());
    solver.solve(system);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // A new right-hand side reuses the retained numeric factorization.
    expected << -2.0, 1.0, 4.0;
    system.getVector() = system.getOperator() * expected;
    EXPECT_TRUE(solver.hasFactorization());
    solver.solve(system);
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // A value-only change keeps the sparsity pattern and reuses analysis, but
    // must rebuild the numeric factorization.
    system.getOperator().coeffRef(0, 0) = 8.0;
    expected << 2.0, -1.0, 4.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_TRUE(solver.hasFactorization());
    solver.solve(system);
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // A structural change requires explicitly clearing symbolic analysis.
    solver.clear(Symbolic);
    system.getOperator().setZero();
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(1, 1) = 3.0;
    system.getOperator().insert(2, 2) = 4.0;
    expected << -1.0, 3.0, 2.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_TRUE(solver.hasFactorization());
    solver.solve(system);
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // Ordering participates in symbolic analysis and therefore invalidates it.
    solver.setOrdering(decltype(solver)::Ordering::Natural);
    EXPECT_FALSE(solver.hasFactorization());
    solver.solve(system);
    EXPECT_TRUE(solver.hasFactorization());
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    solver.clear(Numeric);
    EXPECT_FALSE(solver.hasFactorization());
    solver.solve(system);
    EXPECT_TRUE(solver.hasFactorization());
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  class Rodin_Solver_ParU_AllGeometries : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(Rodin_Solver_ParU_AllGeometries, SolvesVariationalProblemWithCTAD)
  {
    auto mesh = makeMesh(GetParam());

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    RealFunction f = 1.0;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());

    Rodin::Solver::ParU solver(poisson);
    static_assert(
      std::is_same_v<typename FormLanguage::Traits<decltype(solver)>::LinearSystemType,
        typename decltype(poisson)::LinearSystemType>);
    solver.setMaxThreads(1).setOrdering(decltype(solver)::Ordering::AMD).solve();

    const auto& system = poisson.getLinearSystem();
    EXPECT_LT(
      (system.getOperator() * system.getSolution() - system.getVector()).norm(), 1e-11);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, Rodin_Solver_ParU_AllGeometries,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));
}
