/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DenseSolverTest.cpp
 * @brief The solvers that factorize a dense operator.
 *
 * A problem assembles into whichever linear system it is given, so a dense one
 * is named rather than adapted: only the deduction guide of Problem assumes a
 * sparse system, and naming the template arguments bypasses it.
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/HouseholderQR.h"
#include "Rodin/Solver/Info.h"
#include "Rodin/Solver/LDLT.h"
#include "Rodin/Solver/PartialPivLU.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  /// @brief Dense linear system type assembled throughout this file.
  using DenseSystem = Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>;

  /// @brief A problem assembling into a dense operator.
  template <class TrialFunctionType, class TestFunctionType>
  using DenseProblem = Problem<DenseSystem, TrialFunctionType, TestFunctionType>;

  /// @brief A segment mesh with the connectivity a boundary condition needs.
  static Mesh<Context::Local> makeMesh()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {8});
    mesh.getConnectivity().compute(0, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  /// @brief Assembles a Poisson problem into the system the problem carries.
  template <class ProblemType, class TrialFunctionType, class TestFunctionType>
  static void assemblePoisson(
    ProblemType& problem, TrialFunctionType& u, TestFunctionType& v)
  {
    RealFunction f = 1.0;
    problem = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());
    problem.assemble();
  }

  /// @brief Solves and checks that the residual of the assembled system vanishes.
  template <class Solver, class ProblemType>
  static void expectSolvedResidual(Solver& solver, ProblemType& problem)
  {
    solver.solve();

    const auto& system = problem.getLinearSystem();
    EXPECT_GT(system.getOperator().rows(), 0);
    EXPECT_LT(
      (system.getOperator() * system.getSolution() - system.getVector()).norm(), 1e-10);
  }

  /// @brief The dense LDLT solves an assembled problem and reports its status.
  TEST(Rodin_Solver_Dense, LDLTSolvesAnAssembledDenseProblem)
  {
    auto mesh = makeMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    DenseProblem<decltype(u), decltype(v)> poisson(u, v);
    assemblePoisson(poisson, u, v);

    Rodin::Solver::LDLT solver(poisson);
    expectSolvedResidual(solver, poisson);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getInfo().factorization, Rodin::Solver::Factorization::Numeric);
    EXPECT_EQ(solver.getInfo().status, 0);
  }

  /// @brief The dense Householder QR solves an assembled problem.
  TEST(Rodin_Solver_Dense, HouseholderQRSolvesAnAssembledDenseProblem)
  {
    auto mesh = makeMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    DenseProblem<decltype(u), decltype(v)> poisson(u, v);
    assemblePoisson(poisson, u, v);

    // Eigen reports no status for this factorization, so the residual is all
    // there is to assert.
    Rodin::Solver::HouseholderQR solver(poisson);
    expectSolvedResidual(solver, poisson);
  }

  /// @brief The dense partially pivoted LU solves an assembled problem.
  TEST(Rodin_Solver_Dense, PartialPivLUSolvesAnAssembledDenseProblem)
  {
    auto mesh = makeMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    DenseProblem<decltype(u), decltype(v)> poisson(u, v);
    assemblePoisson(poisson, u, v);

    Rodin::Solver::PartialPivLU solver(poisson);
    expectSolvedResidual(solver, poisson);
  }

  /// @brief Dense assembly reaches the same system, and solution, as sparse.
  TEST(Rodin_Solver_Dense, DenseAssemblyMatchesSparseAssembly)
  {
    auto mesh = makeMesh();
    P1 vh(mesh);

    TrialFunction uSparse(vh);
    TestFunction vSparse(vh);
    Problem sparse(uSparse, vSparse);
    assemblePoisson(sparse, uSparse, vSparse);
    Rodin::Solver::SparseLU sparseSolver(sparse);
    sparseSolver.solve();

    TrialFunction uDense(vh);
    TestFunction vDense(vh);
    DenseProblem<decltype(uDense), decltype(vDense)> dense(uDense, vDense);
    assemblePoisson(dense, uDense, vDense);
    Rodin::Solver::LDLT denseSolver(dense);
    denseSolver.solve();

    const auto& sparseSystem = sparse.getLinearSystem();
    const auto& denseSystem = dense.getLinearSystem();
    EXPECT_EQ(denseSystem.getOperator().rows(), sparseSystem.getOperator().rows());
    EXPECT_EQ(denseSystem.getOperator().cols(), sparseSystem.getOperator().cols());
    EXPECT_LT(
      (Math::Matrix<Real>(sparseSystem.getOperator()) - denseSystem.getOperator()).norm(),
      1e-12);
    EXPECT_LT((sparseSystem.getVector() - denseSystem.getVector()).norm(), 1e-12);
    EXPECT_LT((sparseSystem.getSolution() - denseSystem.getSolution()).norm(), 1e-10);
  }
}
