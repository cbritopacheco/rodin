/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
// Exercise the compatibility include; backend tests use the Kokkos header.
#include "Rodin/Solver/CuDSS.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  using Rodin::Solver::Factorization;

  using LinearSystemType =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using CuDSSSolver = Rodin::Solver::CuDSS<LinearSystemType>;
  using Symmetry = CuDSSSolver::Symmetry;

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

  /// @brief A symmetric indefinite system whose solution is (1, 2, 3).
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

  static Math::Vector<Real> oneTwoThree()
  {
    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    return expected;
  }

  static void expectSolution(
    const LinearSystemType& system, const Math::Vector<Real>& expected)
  {
    ASSERT_EQ(system.getSolution().size(), expected.size());
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  static void expectResourceState(const CuDSSSolver& solver, bool symbolic, bool numeric,
    Optional<Factorization> factorization)
  {
    EXPECT_EQ(solver.getResources().symbolic, symbolic);
    EXPECT_EQ(solver.getResources().numeric, numeric);
    EXPECT_EQ(solver.getInfo().factorization, factorization);
  }

  /// @brief A placeholder problem, since the solvers are built from one.
  struct Fixture
  {
    using FESType = decltype(P1(std::declval<const Mesh<Context::Local>&>()));
    using TrialType = decltype(TrialFunction(std::declval<FESType&>()));
    using TestType = decltype(TestFunction(std::declval<FESType&>()));

    Mesh<Context::Local> mesh = makeMesh(Polytope::Type::Segment);
    FESType vh{mesh};
    TrialType u{vh};
    TestType v{vh};
    decltype(Problem(std::declval<TrialType&>(), std::declval<TestType&>())) problem{
      u, v};
  };

  /**
   * @brief The device holds the transpose, so an unsymmetric solve checks
   * that the transpose is undone.
   */
  TEST(Rodin_Solver_CuDSS, SolvesUnsymmetricSystem)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makeUnsymmetricSystem();
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getInfo().status, 0);
    expectSolution(system, oneTwoThree());
    expectResourceState(solver, true, true, Factorization::Numeric);
  }

  TEST(Rodin_Solver_CuDSS, SolvesSymmetricIndefiniteSystem)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makeSymmetricSystem();
    solver.setSymmetric(Symmetry::General).solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());
  }

  TEST(Rodin_Solver_CuDSS, SolvesPositiveDefiniteSystem)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makePositiveDefiniteSystem();
    solver.setSymmetric(Symmetry::PositiveDefinite).solve(system);
    EXPECT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());
  }

  /// @brief New values with the same pattern reuse the symbolic analysis.
  TEST(Rodin_Solver_CuDSS, LifecycleReusesSymbolicAnalysis)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    expectResourceState(solver, false, false, {});

    auto system = makeUnsymmetricSystem();
    solver.factorize(system);
    ASSERT_TRUE(solver.success());
    expectResourceState(solver, true, true, Factorization::Numeric);

    // A second right-hand side reuses the numeric factorization.
    Math::Vector<Real> expected(3);
    expected << -2.0, 0.5, 4.0;
    system.getVector() = system.getOperator() * expected;
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, expected);

    // New values, same pattern.
    system.getOperator().coeffRef(0, 0) = 7.0;
    system.getOperator().coeffRef(2, 1) = -1.0;
    system.getVector() = system.getOperator() * expected;
    solver.factorize(system);
    ASSERT_TRUE(solver.success());
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, expected);

    solver.clear(Factorization::Numeric);
    expectResourceState(solver, true, false, Factorization::Symbolic);
    system.getSolution().setConstant(-99);
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, expected);
    expectResourceState(solver, true, true, Factorization::Numeric);

    solver.clear(Factorization::Symbolic);
    expectResourceState(solver, false, false, {});
    EXPECT_EQ(solver.getResources().values.data, nullptr);

    // A new pattern is accepted once the analysis has been cleared.
    auto other = makePositiveDefiniteSystem();
    solver.solve(other);
    ASSERT_TRUE(solver.success());
    expectSolution(other, oneTwoThree());
  }

  /// @brief A matrix with another pattern is analyzed again.
  TEST(Rodin_Solver_CuDSS, FactorizationAnalyzesNewPattern)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makeUnsymmetricSystem();
    solver.factorize(system);
    ASSERT_TRUE(solver.success());
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());

    // As many entries as the analyzed matrix, but (2, 0) instead of (1, 0).
    LinearSystemType moved;
    moved.getOperator().resize(3, 3);
    moved.getOperator().insert(0, 0) = 4.0;
    moved.getOperator().insert(2, 0) = 2.0;
    moved.getOperator().insert(0, 1) = 1.0;
    moved.getOperator().insert(1, 1) = 5.0;
    moved.getOperator().insert(2, 1) = 3.0;
    moved.getOperator().insert(1, 2) = 1.0;
    moved.getOperator().insert(2, 2) = 6.0;
    moved.getVector().resize(3);
    moved.getVector() << 6.0, 13.0, 26.0;
    solver.factorize(moved);
    ASSERT_TRUE(solver.success());
    expectResourceState(solver, true, true, Factorization::Numeric);
    solver.solve(moved);
    ASSERT_TRUE(solver.success());
    expectSolution(moved, oneTwoThree());

    // Fewer entries.
    auto fewer = makeSymmetricSystem();
    solver.factorize(fewer);
    ASSERT_TRUE(solver.success());
    solver.solve(fewer);
    ASSERT_TRUE(solver.success());
    expectSolution(fewer, oneTwoThree());

    // The first pattern again.
    system.getSolution().setZero();
    solver.factorize(system);
    ASSERT_TRUE(solver.success());
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());
  }

  /// @brief Declaring another symmetry discards the retained factorization.
  TEST(Rodin_Solver_CuDSS, SymmetryChangeReleasesFactorization)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makePositiveDefiniteSystem();
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectResourceState(solver, true, true, Factorization::Numeric);

    solver.setSymmetric(Symmetry::Unsymmetric);
    expectResourceState(solver, true, true, Factorization::Numeric);

    solver.setSymmetric(Symmetry::PositiveDefinite);
    EXPECT_EQ(solver.getSymmetric(), Symmetry::PositiveDefinite);
    expectResourceState(solver, false, false, {});
    system.getSolution().setZero();
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());
  }

  TEST(Rodin_Solver_CuDSS, ReportsUnusableRequests)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makeUnsymmetricSystem();

    // A right-hand side that does not match the matrix is an error, whether
    // or not a factorization is retained.
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

  /// @brief A copy keeps the configuration but owns no device state.
  TEST(Rodin_Solver_CuDSS, CopiesPreserveConfigurationOnly)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makePositiveDefiniteSystem();
    solver.setSymmetric(Symmetry::PositiveDefinite).setRefinementSteps(2).solve(system);
    ASSERT_TRUE(solver.success());

    CuDSSSolver copy(solver);
    EXPECT_EQ(copy.getSymmetric(), Symmetry::PositiveDefinite);
    EXPECT_EQ(copy.getRefinementSteps(), 2u);
    expectResourceState(copy, false, false, {});
    EXPECT_EQ(copy.getResources().handle, nullptr);

    auto again = makePositiveDefiniteSystem();
    copy.solve(again);
    ASSERT_TRUE(copy.success());
    expectSolution(again, oneTwoThree());

    // The original still solves with its own factors.
    system.getSolution().setZero();
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());
  }

  /// @brief Moving configuration must not alias or steal retained factorization state.
  TEST(Rodin_Solver_CuDSS, MovesPreserveConfigurationOnly)
  {
    Fixture f;
    CuDSSSolver solver(f.problem);
    auto system = makePositiveDefiniteSystem();
    solver.setSymmetric(Symmetry::PositiveDefinite).setRefinementSteps(2).solve(system);
    ASSERT_TRUE(solver.success());
    const auto original = solver.getResources().data;

    CuDSSSolver moved(std::move(solver));
    EXPECT_EQ(moved.getSymmetric(), Symmetry::PositiveDefinite);
    EXPECT_EQ(moved.getRefinementSteps(), 2u);
    expectResourceState(moved, false, false, {});
    EXPECT_EQ(moved.getResources().handle, nullptr);
    EXPECT_EQ(solver.getResources().data, original);
    EXPECT_EQ(&moved.getProblem(), &f.problem);

    auto again = makePositiveDefiniteSystem();
    moved.solve(again);
    ASSERT_TRUE(moved.success());
    expectSolution(again, oneTwoThree());
    EXPECT_NE(moved.getResources().data, original);
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    expectSolution(system, oneTwoThree());
  }

  class Rodin_Solver_CuDSS_AllGeometries : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief The solution matches SparseLU on a variational problem.
  TEST_P(Rodin_Solver_CuDSS_AllGeometries, SolvesVariationalProblemWithCTAD)
  {
    auto mesh = makeMesh(GetParam());

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    RealFunction f = 1.0;

    // The advection term makes the operator unsymmetric.
    const auto check = [&](const auto& beta) {
      Problem problem(u, v);
      problem = Integral{Grad{u}, Grad{v}} + Integral{Dot(beta, Grad{u}), v} -
        Integral(f, v) + DirichletBC(u, Zero());

      Rodin::Solver::CuDSS solver(problem);
      static_assert(
        std::is_same_v<typename FormLanguage::Traits<decltype(solver)>::LinearSystemType,
          typename decltype(problem)::LinearSystemType>);
      solver.solve();
      ASSERT_TRUE(solver.success());

      auto& system = problem.getLinearSystem();
      EXPECT_LT(
        (system.getOperator() * system.getSolution() - system.getVector()).norm(), 1e-11);

      const Math::Vector<Real> device = system.getSolution();
      Rodin::Solver::SparseLU reference(problem);
      reference.solve(system);
      ASSERT_TRUE(reference.success());
      EXPECT_LT((device - system.getSolution()).norm(), 1e-10);
    };

    switch (mesh.getSpaceDimension())
    {
      case 1:
      {
        Math::Vector<Real> beta(1);
        beta << 1.0;
        check(VectorFunction(beta));
        break;
      }
      case 2:
        check(VectorFunction{1.0, 0.5});
        break;
      default:
        check(VectorFunction{1.0, 0.5, 0.25});
        break;
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, Rodin_Solver_CuDSS_AllGeometries,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));
}
