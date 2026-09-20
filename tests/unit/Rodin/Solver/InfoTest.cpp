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
  using Rodin::Solver::Factorization;

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

#ifdef RODIN_USE_CHOLMOD
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
