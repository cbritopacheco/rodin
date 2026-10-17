/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InfoTest.cpp
 * @brief The status every factorization-based solver reports.
 *
 * Each of these solvers wraps a library of its own, so the interface is
 * checked on each of them rather than on a shared base.
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/CHOLMOD.h"
#include "Rodin/Solver/Info.h"
#include "Rodin/Solver/SPQR.h"
#include "Rodin/Solver/SimplicialLDLT.h"
#include "Rodin/Solver/SimplicialLLT.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Solver/SparseQR.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  /// @brief Factorization stage reported by a solver.
  using Rodin::Solver::Factorization;

  /// @brief Sparse linear system type solved throughout this file.
  using SparseSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  /// @brief A symmetric positive definite system whose solution is (1, 2, 3).
  static SparseSystem makeSparseSystem()
  {
    SparseSystem system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 4.0;
    system.getOperator().insert(1, 0) = 1.0;
    system.getOperator().insert(0, 1) = 1.0;
    system.getOperator().insert(1, 1) = 5.0;
    system.getOperator().insert(2, 1) = 2.0;
    system.getOperator().insert(1, 2) = 2.0;
    system.getOperator().insert(2, 2) = 6.0;
    system.getOperator().makeCompressed();
    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    system.getVector() = system.getOperator() * expected;
    return system;
  }

  /// @brief A structurally singular system: row and column one are empty.
  static SparseSystem makeSingularSystem()
  {
    SparseSystem system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(2, 2) = 4.0;
    system.getOperator().makeCompressed();
    system.getVector().resize(3);
    system.getVector() << 1.0, 1.0, 1.0;
    // A sentinel the solver must not overwrite with a bogus solution.
    system.getSolution().resize(3);
    system.getSolution() << 7.0, 7.0, 7.0;
    return system;
  }

  template <class Solver, class System>
  static void expectSolved(Solver& solver, System& system)
  {
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getInfo().factorization, Factorization::Numeric);
    EXPECT_EQ(solver.getInfo().status, 0);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  template <class Solver>
  static void expectRefused(Solver& solver, SparseSystem& system)
  {
    solver.solve(system);
    EXPECT_FALSE(solver.success());
    EXPECT_FALSE(solver.getInfo().factorization.has_value());
    EXPECT_NE(solver.getInfo().status, 0);
    // The solve was skipped, so the sentinel survived.
    EXPECT_EQ(system.getSolution()(0), 7.0);
    EXPECT_EQ(system.getSolution()(1), 7.0);
    EXPECT_EQ(system.getSolution()(2), 7.0);
  }

  class Rodin_Solver_Info : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {2});
      }

      Mesh<Context::Local> m_mesh;
  };

  /// @brief A successful factorization is retained and reported.
  TEST_F(Rodin_Solver_Info, SparseQRReportsASuccessfulFactorization)
  {
    P1 vh(m_mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);
    SparseSystem system = makeSparseSystem();

    Rodin::Solver::SparseQR solver(problem);
    expectSolved(solver, system);
  }

  /// @brief A successful factorization is reported, and a failing one refuses to solve.
  TEST_F(Rodin_Solver_Info, SimplicialLLTReportsBothOutcomes)
  {
    P1 vh(m_mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    SparseSystem system = makeSparseSystem();
    Rodin::Solver::SimplicialLLT solver(problem);
    expectSolved(solver, system);

    SparseSystem singular = makeSingularSystem();
    expectRefused(solver, singular);
  }

  /// @brief A successful factorization is reported, and a failing one refuses to solve.
  TEST_F(Rodin_Solver_Info, SimplicialLDLTReportsBothOutcomes)
  {
    P1 vh(m_mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    SparseSystem system = makeSparseSystem();
    Rodin::Solver::SimplicialLDLT solver(problem);
    expectSolved(solver, system);

    SparseSystem singular = makeSingularSystem();
    expectRefused(solver, singular);
  }

  template <class Scalar, template <class> class SolverType>
  void checkOptionalSolverLifecycle()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem problem(u, v);
    using System = Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
    System system;
    system.getOperator().resize(3, 3);
    const Scalar offDiagonal = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, 0.5);
      else
        return Scalar(1);
    }();
    for (Eigen::Index i = 0; i < 3; ++i)
      system.getOperator().insert(i, i) = Scalar(4 + i);
    system.getOperator().insert(0, 1) = offDiagonal;
    system.getOperator().insert(1, 0) = Eigen::numext::conj(offDiagonal);
    system.getOperator().insert(1, 2) = Scalar(-1);
    system.getOperator().insert(2, 1) = Scalar(-1);
    system.getOperator().makeCompressed();
    Math::Vector<Scalar> expected(3);
    expected << Scalar(1), Scalar(-2), Scalar(3);
    system.getVector() = system.getOperator() * expected;
    SolverType<System> solver(problem);
    auto check = [&](auto& candidate) {
      system.getSolution() = Math::Vector<Scalar>::Constant(3, Scalar(7));
      candidate.solve(system);
      EXPECT_TRUE(candidate.success());
      EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
    };
    check(solver);
    auto copied = solver;
    check(copied);
    auto moved = std::move(copied);
    check(moved);
    std::unique_ptr<SolverType<System>> cloned(solver.copy());
    check(*cloned);
  }

#ifdef RODIN_USE_CHOLMOD
  /// @brief Real and complex CHOLMOD solvers retain the problem through copy and move.
  TEST_F(Rodin_Solver_Info, CholmodLifecycle)
  {
    checkOptionalSolverLifecycle<Real, Rodin::Solver::CHOLMOD::SupernodalLLT>();
    checkOptionalSolverLifecycle<Complex, Rodin::Solver::CHOLMOD::SupernodalLLT>();
  }

  /// @brief A successful factorization is retained and reported.
  TEST_F(Rodin_Solver_Info, CholmodSupernodalLLTReportsASuccessfulFactorization)
  {
    P1 vh(m_mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);
    SparseSystem system = makeSparseSystem();

    Rodin::Solver::CHOLMOD::SupernodalLLT solver(problem);
    expectSolved(solver, system);
  }
#endif

#ifdef RODIN_USE_SPQR
  /// @brief Real and complex SPQR solvers retain the problem through copy and move.
  TEST_F(Rodin_Solver_Info, SPQRLifecycle)
  {
    checkOptionalSolverLifecycle<Real, Rodin::Solver::SPQR>();
    checkOptionalSolverLifecycle<Complex, Rodin::Solver::SPQR>();
  }

  /// @brief A successful factorization is retained and reported.
  TEST_F(Rodin_Solver_Info, SPQRReportsASuccessfulFactorization)
  {
    P1 vh(m_mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);
    SparseSystem system = makeSparseSystem();

    Rodin::Solver::SPQR solver(problem);
    expectSolved(solver, system);
  }
#endif
}
