/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/MUMPS.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  using Rodin::Solver::Factorization;

  using LinearSystemType =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

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

  /// @brief An unsymmetric system whose solution is (1, 2, 3).
  static LinearSystemType makeUnsymmetricSystem()
  {
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
    return system;
  }

  /// @brief A symmetric indefinite system, stored with both triangles.
  static LinearSystemType makeSymmetricSystem()
  {
    LinearSystemType system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 1.0;
    system.getOperator().insert(1, 0) = 2.0;
    system.getOperator().insert(0, 1) = 2.0;
    system.getOperator().insert(1, 1) = 1.0;
    system.getOperator().insert(2, 2) = -3.0;
    system.getVector().resize(3);
    system.getVector() << 5.0, 4.0, -9.0;
    return system;
  }

  TEST(Rodin_Solver_MUMPS, SolvesRawLinearSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    Rodin::Solver::MUMPS solver(problem);
    using enum decltype(solver)::Factorization;
    LinearSystemType system = makeUnsymmetricSystem();

    solver.setOrdering(decltype(solver)::Ordering::AMD).factorize(system);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());

    // The system type is process-local, so the solve must be private to this
    // process rather than collective over MPI_COMM_WORLD.
    EXPECT_EQ(solver.getResources().instance.comm_fortran,
      static_cast<MUMPS_INT>(MPI_Comm_c2f(MPI_COMM_SELF)));
    solver.solve(system);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // A new right-hand side reuses the retained numeric factorization.
    expected << -2.0, 1.0, 4.0;
    system.getVector() = system.getOperator() * expected;
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    solver.solve(system);
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // A value-only change keeps the sparsity pattern and reuses analysis, but
    // must rebuild the numeric factorization.
    system.getOperator().coeffRef(0, 0) = 8.0;
    expected << 2.0, -1.0, 4.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    solver.solve(system);
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // A structural change requires explicitly clearing symbolic analysis.
    solver.clear(Symbolic);
    EXPECT_FALSE(solver.getInfo().factorization.has_value());
    system.getOperator().setZero();
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(1, 1) = 3.0;
    system.getOperator().insert(2, 2) = 4.0;
    expected << -1.0, 3.0, 2.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    solver.solve(system);
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    // Ordering participates in symbolic analysis and therefore invalidates it.
    solver.setOrdering(decltype(solver)::Ordering::Automatic);
    EXPECT_NE(solver.getInfo().factorization, Factorization::Numeric);
    solver.solve(system);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);

    solver.clear(Numeric);
    EXPECT_NE(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Symbolic);
    solver.solve(system);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  TEST(Rodin_Solver_MUMPS, SymmetricFactorizationMatchesUnsymmetric)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;

    // The whole matrix is stored; the symmetric factorization must read only
    // its lower triangle and reach the same solution.
    LinearSystemType unsymmetric = makeSymmetricSystem();
    Rodin::Solver::MUMPS reference(problem);
    reference.solve(unsymmetric);
    EXPECT_NEAR((unsymmetric.getSolution() - expected).norm(), 0.0, 1e-12);

    LinearSystemType symmetric = makeSymmetricSystem();
    Rodin::Solver::MUMPS solver(problem);
    solver.setSymmetric(decltype(solver)::Symmetry::General);
    EXPECT_EQ(solver.getSymmetric(), decltype(solver)::Symmetry::General);
    solver.solve(symmetric);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    EXPECT_NEAR((symmetric.getSolution() - expected).norm(), 0.0, 1e-12);

    // Symmetry selects the factorization itself, so changing it releases the
    // retained factorization and the converted structure.
    solver.setSymmetric(decltype(solver)::Symmetry::Unsymmetric);
    EXPECT_NE(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_FALSE(solver.getInfo().factorization.has_value());
    solver.solve(symmetric);
    EXPECT_NEAR((symmetric.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  TEST(Rodin_Solver_MUMPS, SolvesPositiveDefiniteSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    LinearSystemType system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 4.0;
    system.getOperator().insert(1, 0) = 1.0;
    system.getOperator().insert(0, 1) = 1.0;
    system.getOperator().insert(1, 1) = 5.0;
    system.getOperator().insert(2, 1) = 2.0;
    system.getOperator().insert(1, 2) = 2.0;
    system.getOperator().insert(2, 2) = 6.0;
    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    system.getVector() = system.getOperator() * expected;

    Rodin::Solver::MUMPS solver(problem);
    solver.setSymmetric(decltype(solver)::Symmetry::PositiveDefinite).solve(system);
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_TRUE(solver.success());
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  TEST(Rodin_Solver_MUMPS, ReportsUnusableRequests)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    Rodin::Solver::MUMPS solver(problem);
    LinearSystemType system = makeUnsymmetricSystem();

    // A right-hand side that does not match the matrix is an error, whether or
    // not a factorization is retained.
    system.getVector().resize(2);
    EXPECT_ANY_THROW(solver.solve(system));
    system = makeUnsymmetricSystem();
    solver.factorize(system);
    system.getVector().resize(2);
    EXPECT_ANY_THROW(solver.solve(system));

    // A rectangular matrix cannot be factorized.
    LinearSystemType rectangular;
    rectangular.getOperator().resize(2, 3);
    rectangular.getVector().resize(2);
    rectangular.getVector().setZero();
    EXPECT_ANY_THROW(solver.factorize(rectangular));
  }

  class Rodin_Solver_MUMPS_AllGeometries : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(Rodin_Solver_MUMPS_AllGeometries, SolvesVariationalProblemWithCTAD)
  {
    auto mesh = makeMesh(GetParam());

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    RealFunction f = 1.0;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());

    Rodin::Solver::MUMPS solver(poisson);
    static_assert(
      std::is_same_v<typename FormLanguage::Traits<decltype(solver)>::LinearSystemType,
        typename decltype(poisson)::LinearSystemType>);
    solver.setOrdering(decltype(solver)::Ordering::AMD).solve();

    const auto& system = poisson.getLinearSystem();
    EXPECT_LT(
      (system.getOperator() * system.getSolution() - system.getVector()).norm(), 1e-11);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, Rodin_Solver_MUMPS_AllGeometries,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));
}
