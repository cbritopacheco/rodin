/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <memory>
#include <cmath>
#include <string>
#include "Rodin/Variational/Problem.h"
#include "Rodin/Adaptation/SWIFT/LinearSolver.h"

namespace Rodin::Tests::Unit
{
  namespace SWIFT = Adaptation::SWIFT;
  using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  class SWIFTLinearSolverTest : public ::testing::Test
  {
    protected:
      class Adapter final : public Variational::ProblemBase<System>
      {
        public:
          Adapter& operator=(const ProblemBodyType&) override
          {
            return *this;
          }
          Adapter& assemble() override
          {
            return *this;
          }
          void solve(Solver::LinearSolverBase<System>& solver) override
          {
            solver.solve(m_system);
          }
          System& getLinearSystem() override
          {
            return m_system;
          }
          const System& getLinearSystem() const override
          {
            return m_system;
          }
          Adapter* copy() const noexcept override
          {
            return new Adapter(*this);
          }

        private:
          System m_system;
      };

      void SetUp() override
      {
        Math::Matrix<Real> matrix(2, 2);
        matrix << 4, 1, 1, 3;
        auto& system = problem.getLinearSystem();
        system.getOperator() = matrix.sparseView();
        system.getOperator().makeCompressed();
        system.getVector().resize(2);
        system.getVector() << Real(5.5), Real(6.75);
        system.getSolution() = Math::Vector<Real>::Zero(2);
        centering.resize(2, 1);
        centering << Real(0.5), Real(0.25);
        parameters.linear.solver = SWIFT::Parameters::LinearSolver::CG;
        parameters.convergence.iterations.linear = 0;
        parameters.convergence.tolerance.linearRelative = tolerance;
      }

      void expectCenteredSolution(const SWIFT::LinearSolver<System>& solver)
      {
        const auto& system = problem.getLinearSystem();
        ASSERT_TRUE(solver.success());
        ASSERT_EQ(system.getSolution().size(), 2);
        EXPECT_NEAR(system.getSolution()(0), 1, tolerance);
        EXPECT_NEAR(system.getSolution()(1), 2, tolerance);
        Math::Vector<Real> residual = system.getOperator() * system.getSolution() -
          centering * (centering.transpose() * system.getSolution()) - system.getVector();
        EXPECT_LE(residual.norm(), tolerance * system.getVector().norm());
        EXPECT_LE(solver.getError(), tolerance);
      }

      /// Roundoff allowance for exact two-dimensional algebraic solutions.
      static constexpr Real tolerance = 1e-12;
      Adapter problem;
      Math::Matrix<Real> centering;
      SWIFT::Parameters parameters;
  };

  /// @brief A copied solver retains its problem, centered operator and selected backend.
  TEST_F(SWIFTLinearSolverTest, CopiesOwnConfigurationAndRemainBoundToProblem)
  {
    for (auto backend :
      {SWIFT::Parameters::LinearSolver::CG, SWIFT::Parameters::LinearSolver::SparseLU})
    {
      SCOPED_TRACE(static_cast<int>(backend));
      auto& system = problem.getLinearSystem();
      system.getSolution().setZero();
      parameters.linear.solver = backend;
      SWIFT::LinearSolver<System> original(problem);
      original.setParameters(parameters).setCentering(centering);
      original.solve();
      expectCenteredSolution(original);
      std::unique_ptr<SWIFT::LinearSolver<System>> cloned(original.copy());
      EXPECT_EQ(&cloned->getProblem(), &problem);
      auto altered = parameters;
      altered.linear.solver = SWIFT::Parameters::LinearSolver::CG;
      altered.convergence.iterations.linear = 1;
      original.setParameters(altered).setCentering(Math::Matrix<Real>());
      system.getSolution().setZero();
      cloned->solve();
      expectCenteredSolution(*cloned);
      if (backend == SWIFT::Parameters::LinearSolver::CG)
        EXPECT_EQ(cloned->getIterations(), 2u);
      else
        EXPECT_EQ(cloned->getIterations(), 0u);
      original.solve();
      EXPECT_FALSE(original.success());
      system.getSolution().setZero();
      cloned->solve();
      expectCenteredSolution(*cloned);
    }
  }

  /// @brief Zero forcing overwrites a stale guess with the exact zero solution.
  TEST_F(SWIFTLinearSolverTest, ZeroLoadClearsGuessWithNoIterations)
  {
    auto& system = problem.getLinearSystem();
    const auto matrix = system.getOperator();
    system.getVector().setZero();
    system.getSolution().setConstant(7);
    SWIFT::LinearSolver<System> solver(problem);
    solver.setParameters(parameters).setCentering(centering).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterations(), 0u);
    EXPECT_EQ(solver.getError(), 0);
    EXPECT_EQ(system.getSolution().norm(), 0);
    EXPECT_EQ((system.getOperator() - matrix).norm(), 0);
  }

  /// @brief An exact warm start terminates before the first correction or division.
  TEST_F(SWIFTLinearSolverTest, ExactWarmStartNeedsNoCorrection)
  {
    auto& system = problem.getLinearSystem();
    system.getSolution() << 1, 2;
    SWIFT::LinearSolver<System> solver(problem);
    solver.setParameters(parameters).setCentering(centering).solve();
    expectCenteredSolution(solver);
    EXPECT_EQ(solver.getIterations(), 0u);
  }

  /// @brief Nonpositive curvature fails finitely, then the same solver accepts an SPD system.
  TEST_F(SWIFTLinearSolverTest, NegativeCurvatureReportsFailureAndRecovers)
  {
    auto& system = problem.getLinearSystem();
    const auto positive = system.getOperator();
    const auto load = system.getVector();
    Math::Matrix<Real> indefinite(2, 2);
    indefinite << 1, 2, 2, 1;
    system.getOperator() = indefinite.sparseView();
    system.getVector() << 1, -1;
    parameters.trace = true;
    SWIFT::LinearSolver<System> solver(problem);
    solver.setParameters(parameters);
    testing::internal::CaptureStdout();
    solver.solve();
    const auto diagnostic = testing::internal::GetCapturedStdout();
    EXPECT_NE(diagnostic.find("cg failure"), std::string::npos);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterations(), 0u);
    EXPECT_EQ(solver.getError(), 1);
    EXPECT_EQ(system.getSolution().norm(), 0);
    EXPECT_EQ((system.getOperator() - indefinite.sparseView()).norm(), 0);
    system.getOperator() = positive;
    system.getVector() = load;
    solver.setCentering(centering).solve();
    expectCenteredSolution(solver);
  }

  /// @brief Unsupported targeted assembly reports member context and preserves all system data.
  TEST_F(SWIFTLinearSolverTest, UnsupportedAssemblyTargetPreservesBoundSystem)
  {
    Variational::ProblemBase<System>& base = problem;
    auto& system = base.getLinearSystem();
    const auto matrix = system.getOperator();
    const auto load = system.getVector();
    const auto guess = system.getSolution();
    for (auto target :
      {Variational::AssemblyTarget::LHS, Variational::AssemblyTarget::RHS})
    {
      try
      {
        base.assemble(target);
        ADD_FAILURE() << "Adapter accepted unsupported targeted assembly";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("member function assemble"), std::string::npos);
        EXPECT_NE(
          message.find("Targeted assembly is not implemented"), std::string::npos);
      }
      EXPECT_EQ((system.getOperator() - matrix).norm(), 0);
      EXPECT_EQ((system.getVector() - load).norm(), 0);
      EXPECT_EQ((system.getSolution() - guess).norm(), 0);
    }
    SWIFT::LinearSolver<System> solver(problem);
    solver.setParameters(parameters).setCentering(centering).solve();
    expectCenteredSolution(solver);
  }

  /// @brief Singular direct factorization preserves data and the same adapter can solve again.
  TEST_F(SWIFTLinearSolverTest, SingularDirectFactorizationReportsFailureAndRecovers)
  {
    auto& system = problem.getLinearSystem();
    const auto positive = system.getOperator();
    const auto load = system.getVector();
    Math::Matrix<Real> singular(2, 2);
    singular << 1, 2, 2, 4;
    system.getOperator() = singular.sparseView();
    system.getOperator().makeCompressed();
    system.getSolution().setConstant(17);
    const auto guess = system.getSolution();
    parameters.linear.solver = SWIFT::Parameters::LinearSolver::SparseLU;
    parameters.trace = true;
    SWIFT::LinearSolver<System> solver(problem);
    solver.setParameters(parameters);
    testing::internal::CaptureStdout();
    solver.solve();
    const auto diagnostic = testing::internal::GetCapturedStdout();
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterations(), 0u);
    EXPECT_TRUE(std::isinf(solver.getError()));
    EXPECT_GT(solver.getError(), 0);
    EXPECT_NE(diagnostic.find("metric LU: ok=0"), std::string::npos);
    EXPECT_EQ((system.getSolution() - guess).norm(), 0);
    EXPECT_EQ((system.getVector() - load).norm(), 0);
    EXPECT_EQ((system.getOperator() - singular.sparseView()).norm(), 0);
    system.getOperator() = positive;
    solver.setCentering(centering).solve();
    expectCenteredSolution(solver);
    EXPECT_EQ(solver.getIterations(), 0u);
  }

}
