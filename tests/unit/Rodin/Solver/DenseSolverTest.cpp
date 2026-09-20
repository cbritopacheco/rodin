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
 * Assembly produces sparse systems, so no problem in the tree carries a dense
 * one. These solvers are exercised through a problem that holds a dense system
 * and does nothing else, which is all they read from it.
 */
#include <gtest/gtest.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/HouseholderQR.h"
#include "Rodin/Solver/Info.h"
#include "Rodin/Solver/LDLT.h"
#include "Rodin/Solver/PartialPivLU.h"
#include "Rodin/Variational/Problem.h"

using namespace Rodin;

namespace Rodin::Tests::Unit::Solver
{
  using DenseSystem = Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>;

  /// @brief A problem that carries a dense system and assembles nothing.
  class DenseProblem final : public Variational::ProblemBase<DenseSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      DenseProblem& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Rodin::Solver::LinearSolverBase<DenseSystem>& solver) override
      {
        solver.solve(m_system);
      }

      DenseProblem& assemble() override
      {
        return *this;
      }

      DenseSystem& getLinearSystem() override
      {
        return m_system;
      }

      const DenseSystem& getLinearSystem() const override
      {
        return m_system;
      }

      DenseProblem* copy() const noexcept override
      {
        return new DenseProblem(*this);
      }

    private:
      DenseSystem m_system;
  };

  /// @brief A symmetric positive definite system whose solution is (1, 2, 3).
  static void fill(DenseSystem& system)
  {
    system.getOperator().resize(3, 3);
    system.getOperator() << 4.0, 1.0, 0.0, 1.0, 5.0, 2.0, 0.0, 2.0, 6.0;
    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    system.getVector() = system.getOperator() * expected;
  }

  template <class Solver>
  static void expectSolved(Solver& solver, DenseSystem& system)
  {
    solver.solve(system);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  TEST(Rodin_Solver_Dense, LDLTSolvesAndReportsItsStatus)
  {
    DenseProblem problem;
    fill(problem.getLinearSystem());

    Rodin::Solver::LDLT solver(problem);
    expectSolved(solver, problem.getLinearSystem());
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getInfo().factorization, Rodin::Solver::Factorization::Numeric);
    EXPECT_EQ(solver.getInfo().status, 0);
  }

  TEST(Rodin_Solver_Dense, HouseholderQRSolves)
  {
    DenseProblem problem;
    fill(problem.getLinearSystem());

    // Eigen reports no status for this factorization, so there is nothing to
    // assert beyond the solution itself.
    Rodin::Solver::HouseholderQR solver(problem);
    expectSolved(solver, problem.getLinearSystem());
  }

  TEST(Rodin_Solver_Dense, PartialPivLUSolves)
  {
    DenseProblem problem;
    fill(problem.getLinearSystem());

    Rodin::Solver::PartialPivLU solver(problem);
    expectSolved(solver, problem.getLinearSystem());
  }

  /// @brief Solving through the problem reaches the same solution.
  TEST(Rodin_Solver_Dense, SolvesThroughTheProblem)
  {
    DenseProblem problem;
    fill(problem.getLinearSystem());

    Rodin::Solver::LDLT solver(problem);
    solver.solve();

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((problem.getLinearSystem().getSolution() - expected).norm(), 0.0, 1e-12);
  }
}
