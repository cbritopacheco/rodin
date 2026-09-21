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
  using MUMPSSolver = Rodin::Solver::MUMPS<LinearSystemType>;

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

  /// @brief A symmetric positive definite system whose solution is (1, 2, 3).
  static LinearSystemType makePositiveDefiniteSystem()
  {
    LinearSystemType system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 4.0;
    system.getOperator().insert(1, 0) = 1.0;
    system.getOperator().insert(0, 1) = 1.0;
    system.getOperator().insert(1, 1) = 5.0;
    system.getOperator().insert(2, 1) = 2.0;
    system.getOperator().insert(1, 2) = 2.0;
    system.getOperator().insert(2, 2) = 6.0;
    system.getVector().resize(3);
    system.getVector() << 6.0, 17.0, 22.0;
    return system;
  }

  static void expectSolution(const LinearSystemType& system, const Math::Vector<Real>& expected)
  {
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  static void expectResourceState(const MUMPSSolver& solver, bool initialized, bool symbolic,
    bool numeric, Optional<Factorization> factorization)
  {
    EXPECT_EQ(solver.getResources().initialized, initialized);
    EXPECT_EQ(solver.getResources().symbolic, symbolic);
    EXPECT_EQ(solver.getResources().numeric, numeric);
    EXPECT_EQ(solver.getInfo().factorization, factorization);
  }

  TEST(Rodin_Solver_MUMPS, UnsymmetricLifecycleTracksRetainedStages)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    MUMPSSolver solver(problem);
    using enum decltype(solver)::Factorization;
    LinearSystemType system = makeUnsymmetricSystem();
    Math::Vector<Real> expected(3);

    expectResourceState(solver, false, false, false, {});

    solver.setOrdering(decltype(solver)::Ordering::AMD).factorize(system);
    expectResourceState(solver, true, true, true, Factorization::Numeric);
    EXPECT_TRUE(solver.success());

    // The system type is process-local, so the solve must be private to this
    // process rather than collective over MPI_COMM_WORLD.
    EXPECT_EQ(solver.getResources().instance.comm_fortran,
      static_cast<MUMPS_INT>(MPI_Comm_c2f(MPI_COMM_SELF)));
    solver.solve(system);

    expected << 1.0, 2.0, 3.0;
    expectSolution(system, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    // A new right-hand side reuses the retained numeric factorization.
    expected << -2.0, 1.0, 4.0;
    system.getVector() = system.getOperator() * expected;
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    // A value-only change keeps the sparsity pattern and reuses analysis, but
    // must rebuild the numeric factorization.
    system.getOperator().coeffRef(0, 0) = 8.0;
    expected << 2.0, -1.0, 4.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_TRUE(solver.success());
    expectResourceState(solver, true, true, true, Factorization::Numeric);
    solver.solve(system);
    expectSolution(system, expected);

    solver.clear(Numeric);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);
    EXPECT_EQ(solver.getResources().instance.comm_fortran,
      static_cast<MUMPS_INT>(MPI_Comm_c2f(MPI_COMM_SELF)));

    solver.clear(Numeric);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);

    solver.solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    // A structural change requires explicitly clearing symbolic analysis.
    solver.clear(Symbolic);
    expectResourceState(solver, false, false, false, {});
    system.getOperator().setZero();
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(1, 1) = 3.0;
    system.getOperator().insert(2, 2) = 4.0;
    expected << -1.0, 3.0, 2.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_TRUE(solver.success());
    expectResourceState(solver, true, true, true, Factorization::Numeric);
    solver.solve(system);
    expectSolution(system, expected);

    solver.clear(Numeric);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);
    solver.clear(Symbolic);
    expectResourceState(solver, false, false, false, {});
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    // Ordering participates in symbolic analysis and therefore invalidates it.
    solver.setOrdering(decltype(solver)::Ordering::Automatic);
    expectResourceState(solver, false, false, false, {});
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);
  }

  TEST(Rodin_Solver_MUMPS, SymmetricLifecycleAndSymmetryChanges)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);
    using enum Factorization;

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;

    // The whole matrix is stored; the symmetric factorization must read only
    // its lower triangle and reach the same solution.
    LinearSystemType unsymmetric = makeSymmetricSystem();
    Rodin::Solver::MUMPS reference(problem);
    reference.solve(unsymmetric);
    EXPECT_NEAR((unsymmetric.getSolution() - expected).norm(), 0.0, 1e-12);

    LinearSystemType symmetric = makeSymmetricSystem();
    MUMPSSolver solver(problem);
    solver.setSymmetric(decltype(solver)::Symmetry::General);
    EXPECT_EQ(solver.getSymmetric(), decltype(solver)::Symmetry::General);
    expectResourceState(solver, false, false, false, {});
    solver.solve(symmetric);
    EXPECT_TRUE(solver.success());
    expectSolution(symmetric, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    solver.clear(Numeric);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);
    solver.clear(Numeric);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);
    solver.solve(symmetric);
    EXPECT_TRUE(solver.success());
    expectSolution(symmetric, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    solver.clear(Numeric);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);
    solver.clear(Symbolic);
    expectResourceState(solver, false, false, false, {});
    solver.solve(symmetric);
    EXPECT_TRUE(solver.success());
    expectSolution(symmetric, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);

    // Symmetry selects the factorization itself, so changing it releases the
    // retained factorization and the converted structure.
    solver.setSymmetric(decltype(solver)::Symmetry::Unsymmetric);
    expectResourceState(solver, false, false, false, {});
    solver.solve(symmetric);
    EXPECT_TRUE(solver.success());
    expectSolution(symmetric, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);
  }

  TEST(Rodin_Solver_MUMPS, SolvesPositiveDefiniteSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    LinearSystemType system = makePositiveDefiniteSystem();
    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;

    MUMPSSolver solver(problem);
    solver.setSymmetric(decltype(solver)::Symmetry::PositiveDefinite).solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, expected);
    expectResourceState(solver, true, true, true, Factorization::Numeric);
  }

  TEST(Rodin_Solver_MUMPS, FailedRefactorizationNeverReusesInvalidNumericFactors)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    MUMPSSolver solver(problem);
    using enum decltype(solver)::Factorization;
    LinearSystemType system = makePositiveDefiniteSystem();
    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;

    solver.setSymmetric(decltype(solver)::Symmetry::PositiveDefinite).factorize(system);
    ASSERT_TRUE(solver.success());
    expectResourceState(solver, true, true, true, Factorization::Numeric);
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, expected);

    // Same sparsity pattern, but now singular: symbolic analysis remains valid,
    // whereas the new numeric factorization must fail and must not leave the
    // old factors available.
    system.getOperator().coeffRef(0, 0) = 1.0;
    system.getOperator().coeffRef(1, 0) = 1.0;
    system.getOperator().coeffRef(0, 1) = 1.0;
    system.getOperator().coeffRef(1, 1) = 2.0;
    system.getOperator().coeffRef(2, 1) = 1.0;
    system.getOperator().coeffRef(1, 2) = 1.0;
    system.getOperator().coeffRef(2, 2) = 1.0;
    expected << 2.0, -1.0, 3.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    EXPECT_FALSE(solver.success());
    EXPECT_NE(solver.getInfo().status, 0);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);

    system.getSolution().resize(3);
    system.getSolution() << 7.0, 7.0, 7.0;
    solver.solve(system);
    EXPECT_FALSE(solver.success());
    EXPECT_NE(solver.getInfo().status, 0);
    expectResourceState(solver, true, true, false, Factorization::Symbolic);
    EXPECT_EQ(system.getSolution()(0), 7.0);
    EXPECT_EQ(system.getSolution()(1), 7.0);
    EXPECT_EQ(system.getSolution()(2), 7.0);
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

  /// @brief A problem of several trial and test functions deduces the solver.
  TEST(Rodin_Solver_MUMPS, DeducesProblemOfSeveralFields)
  {
    auto mesh = makeMesh(Polytope::Type::Triangle);

    P1 vh(mesh);
    P1 qh(mesh);
    TrialFunction u(vh);
    TrialFunction p(qh);
    TestFunction v(vh);
    TestFunction q(qh);
    RealFunction f = 1.0;

    Problem mixed(u, p, v, q);
    mixed = Integral(Grad(u), Grad(v)) + Integral(p, q) - Integral(f, v) -
      Integral(f, q) + DirichletBC(u, Zero());

    Rodin::Solver::MUMPS solver(mixed);
    static_assert(
      std::is_same_v<typename FormLanguage::Traits<decltype(solver)>::LinearSystemType,
        typename decltype(mixed)::LinearSystemType>);
    solver.solve();
    EXPECT_TRUE(solver.success());

    const auto& system = mixed.getLinearSystem();
    EXPECT_LT(
      (system.getOperator() * system.getSolution() - system.getVector()).norm(), 1e-11);
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
