/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Variational/Problem.h"

using namespace Rodin;

namespace Rodin::Tests::Unit::NewtonSolverTestHelpers
{
  using DenseLinearSystem = Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>;

  class DenseLinearSolver final : public Solver::LinearSolverBase<DenseLinearSystem>
  {
    public:
      using Parent = Solver::LinearSolverBase<DenseLinearSystem>;
      using ProblemBaseType = typename Parent::ProblemBaseType;

      using Parent::solve;

      explicit DenseLinearSolver(ProblemBaseType& pb, bool nonfinite = false)
        : Parent(pb),
          m_nonfinite(nonfinite)
      {}

      void solve(DenseLinearSystem& system) override
      {
        if (m_nonfinite)
        {
          system.getSolution() =
            Math::Vector<Real>::Constant(1, std::numeric_limits<Real>::quiet_NaN());
        }
        else
        {
          system.getSolution() =
            system.getOperator().fullPivLu().solve(system.getVector());
        }
      }

      DenseLinearSolver* copy() const noexcept override
      {
        return new DenseLinearSolver(*this);
      }

    private:
      /// @brief Whether to inject a nonfinite correction for failure-path checks.
      bool m_nonfinite;
  };

  class ScalarNonlinearProblem final : public Variational::ProblemBase<DenseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      explicit ScalarNonlinearProblem(Math::Vector<Real>& u)
        : m_solution(u)
      {}

      Parent& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Solver::LinearSolverBase<DenseLinearSystem>& solver) override
      {
        assemble();
        solver.solve(m_system);
      }

      ScalarNonlinearProblem& assemble() override
      {
        m_system.getOperator().resize(1, 1);
        m_system.getVector().resize(1);
        m_system.getSolution().resize(1);

        const Real u = m_solution.get()(0);
        m_system.getOperator()(0, 0) = 2.0 * u;
        m_system.getVector()(0) = 2.0 - u * u;
        return *this;
      }

      DenseLinearSystem& getLinearSystem() override
      {
        return m_system;
      }

      const DenseLinearSystem& getLinearSystem() const override
      {
        return m_system;
      }

      ScalarNonlinearProblem* copy() const noexcept override
      {
        return new ScalarNonlinearProblem(*this);
      }

    private:
      std::reference_wrapper<Math::Vector<Real>> m_solution;
      DenseLinearSystem m_system;
  };

  class FailingAssembleProblem final : public Variational::ProblemBase<DenseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      Parent& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Solver::LinearSolverBase<DenseLinearSystem>&) override
      {}

      FailingAssembleProblem& assemble() override
      {
        throw std::runtime_error("Assembly failed");
      }

      DenseLinearSystem& getLinearSystem() override
      {
        return m_system;
      }

      const DenseLinearSystem& getLinearSystem() const override
      {
        return m_system;
      }

      FailingAssembleProblem* copy() const noexcept override
      {
        return new FailingAssembleProblem(*this);
      }

    private:
      DenseLinearSystem m_system;
  };
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<Tests::Unit::NewtonSolverTestHelpers::DenseLinearSolver>
  {
    using LinearSystemType = Tests::Unit::NewtonSolverTestHelpers::DenseLinearSystem;
  };
}

using namespace Rodin::Tests::Unit::NewtonSolverTestHelpers;

/// @brief Verifies solves scalar problem using problem assembly for newton solver test by checking tolerance-based numerical results, solver behavior.
TEST(NewtonSolverTest, SolvesScalarProblemUsingProblemAssembly)
{
  Math::Vector<Real> u(1);
  u << 1.5;

  ScalarNonlinearProblem pb(u);
  DenseLinearSolver linearSolver(pb);
  Solver::NewtonSolver newton(linearSolver);
  newton.setMaxIterations(30)
    .setAbsoluteTolerance(1e-12)
    .setRelativeTolerance(1e-12);

  newton.solve(u);

  EXPECT_NEAR(u(0), std::sqrt(2.0), 1e-10);
}

/// @brief Verifies scalar residual shows quadratic convergence rate for newton solver test by checking tolerance-based numerical results, true predicates, solver behavior.
TEST(NewtonSolverTest, ScalarResidualShowsQuadraticConvergenceRate)
{
  Math::Vector<Real> u(1);
  u << 1.45;

  ScalarNonlinearProblem pb(u);
  DenseLinearSolver linearSolver(pb);
  Solver::NewtonSolver newton(linearSolver);
  std::vector<Real> residuals;

  newton.setMaxIterations(20)
    .setAbsoluteTolerance(1e-15)
    .setRelativeTolerance(1e-15)
    .setStepTolerance(0.0)
    .setMonitor([&](const auto& report) { residuals.push_back(report.finalResidual); });

  newton.solve(u);

  ASSERT_GE(residuals.size(), 4u);
  // Skip the first transition and check the asymptotic quadratic regime.
  for (size_t i = 1; i + 1 < residuals.size(); ++i)
  {
    const Real r0 = residuals[i];
    const Real r1 = residuals[i + 1];
    if (!(r0 > 1e-14))
      continue;
    const Real q = r1 / (r0 * r0);
    EXPECT_TRUE(std::isfinite(q));
    EXPECT_GT(q, 0.0);
    EXPECT_LT(q, 0.5);
  }

  EXPECT_NEAR(u(0), std::sqrt(2.0), 1e-12);
}

/// @brief Verifies single template parameter deduces linear system for newton solver test by checking solver behavior.
TEST(NewtonSolverTest, SingleTemplateParameterDeducesLinearSystem)
{
  using SolverType = Solver::NewtonSolver<DenseLinearSolver>;
  static_assert(std::is_same_v<SolverType::LinearSystemType, DenseLinearSystem>);
  SUCCEED();
}

/// @brief Verifies propagates assembly failure for newton solver test by checking exception behavior, solver behavior.
TEST(NewtonSolverTest, PropagatesAssemblyFailure)
{
  FailingAssembleProblem pb;
  DenseLinearSolver linearSolver(pb);
  Solver::NewtonSolver newton(linearSolver);

  Math::Vector<Real> u(1);
  u << 1.0;
  EXPECT_THROW(newton.solve(u), std::runtime_error);
}

/// @brief Invalid solver controls are rejected without overwriting accepted settings.
TEST(NewtonSolverTest, RejectsInvalidControls)
{
  Math::Vector<Real> state(1);
  state << 1;
  ScalarNonlinearProblem problem(state);
  DenseLinearSolver linear(problem);
  Solver::NewtonSolver newton(linear);
  for (Real invalid :
    {-1.0, std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW(newton.setAbsoluteTolerance(invalid), Alert::Exception);
    EXPECT_THROW(newton.setRelativeTolerance(invalid), Alert::Exception);
    EXPECT_THROW(newton.setStepTolerance(invalid), Alert::Exception);
    EXPECT_THROW(newton.setDampingFactor(invalid), Alert::Exception);
  }
  EXPECT_THROW(newton.setDampingFactor(0), Alert::Exception);
  EXPECT_THROW(newton.setMaxIterations(0), Alert::Exception);
  EXPECT_EQ(newton.getAbsoluteTolerance(), 1e-12);
  EXPECT_EQ(newton.getRelativeTolerance(), 1e-8);
  EXPECT_EQ(newton.getStepTolerance(), 0);
  EXPECT_EQ(newton.getDampingFactor(), 1);
  EXPECT_EQ(newton.getMaxIterations(), 100);
  newton.setAbsoluteTolerance(0).setRelativeTolerance(0).setStepTolerance(0);
  EXPECT_EQ(newton.getAbsoluteTolerance(), 0);
  EXPECT_EQ(newton.getRelativeTolerance(), 0);
}

/// @brief A single damped Newton correction has an independently known state and norm.
TEST(NewtonSolverTest, DampedStepAndIterationLimit)
{
  Math::Vector<Real> state(1);
  state << 1;
  ScalarNonlinearProblem problem(state);
  DenseLinearSolver linear(problem);
  Solver::NewtonSolver newton(linear);
  newton.setDampingFactor(0.5).setMaxIterations(1);
  newton.solve(state);
  EXPECT_NEAR(state[0], 1.25, 1e-14);
  const auto& report = newton.getReport();
  EXPECT_FALSE(newton.converged());
  EXPECT_EQ(report.reason, decltype(newton)::ConvergedReason::MaxIterations);
  EXPECT_EQ(report.iterations, 1);
  EXPECT_NEAR(report.initialResidual, 1, 1e-14);
  EXPECT_NEAR(report.finalStepNorm, 0.25, 1e-14);
  EXPECT_EQ(report.dampingFactor, 0.5);
}

/// @brief Relative termination uses the initial residual rather than an absolute threshold.
TEST(NewtonSolverTest, RelativeResidualTermination)
{
  Math::Vector<Real> state(1);
  state << 1;
  ScalarNonlinearProblem problem(state);
  DenseLinearSolver linear(problem);
  Solver::NewtonSolver newton(linear);
  newton.setAbsoluteTolerance(0).setRelativeTolerance(0.3);
  newton.solve(state);
  EXPECT_NEAR(state[0], 1.5, 1e-14);
  EXPECT_TRUE(newton.converged());
  EXPECT_EQ(
    newton.getReport().reason, decltype(newton)::ConvergedReason::RelativeTolerance);
  EXPECT_NEAR(newton.getReport().finalResidual, 0.25, 1e-14);
  EXPECT_EQ(newton.getReport().iterations, 1);
}

/// @brief Step termination applies the accepted increment before returning.
TEST(NewtonSolverTest, StepToleranceAppliesCorrection)
{
  Math::Vector<Real> state(1);
  state << 1;
  ScalarNonlinearProblem problem(state);
  DenseLinearSolver linear(problem);
  Solver::NewtonSolver newton(linear);
  newton.setStepTolerance(0.6);
  newton.solve(state);
  EXPECT_NEAR(state[0], 1.5, 1e-14);
  EXPECT_TRUE(newton.converged());
  EXPECT_EQ(newton.getReport().reason, decltype(newton)::ConvergedReason::StepTolerance);
  EXPECT_NEAR(newton.getReport().finalStepNorm, 0.5, 1e-14);
}

/// @brief Nonfinite residuals or linear corrections leave the initial state unchanged.
TEST(NewtonSolverTest, NonfiniteResidualAndCorrection)
{
  for (bool residualFailure : {false, true})
  {
    Math::Vector<Real> state(1);
    state << (residualFailure ? std::numeric_limits<Real>::infinity() : 1);
    ScalarNonlinearProblem problem(state);
    DenseLinearSolver linear(problem, !residualFailure);
    Solver::NewtonSolver newton(linear);
    size_t notifications = 0;
    newton.setMonitor([&](const auto&) { ++notifications; });
    newton.solve(state);
    using Reason = decltype(newton)::ConvergedReason;
    EXPECT_EQ(newton.getReport().reason,
      residualFailure ? Reason::ResidualNormIsNotFinite : Reason::StepNormIsNotFinite);
    EXPECT_FALSE(newton.converged());
    EXPECT_EQ(notifications, 1);
    if (residualFailure)
      EXPECT_TRUE(std::isinf(state[0]));
    else
      EXPECT_EQ(state[0], 1);
  }
}

/// @brief Custom policies control acceptance, user termination, nonfinite steps and state updates.
TEST(NewtonSolverTest, CustomStepPolicyOutcomes)
{
  for (size_t mode = 0; mode < 4; ++mode)
  {
    Math::Vector<Real> state(1);
    state << 1;
    ScalarNonlinearProblem problem(state);
    DenseLinearSolver linear(problem);
    Solver::NewtonSolver newton(linear);
    newton.setRelativeTolerance(0).setAbsoluteTolerance(1e-14);
    using Newton = decltype(newton);
    size_t calls = 0;
    newton.setStepPolicy([&](auto& x, auto& system, auto&) {
      ++calls;
      typename Newton::StepResult result;
      result.stepNorm = system.getSolution().norm();
      if (mode == 0)
        result.accepted = false;
      else if (mode == 1)
        result.converged = true;
      else if (mode == 2)
        result.stepNorm = std::numeric_limits<Real>::quiet_NaN();
      else
        x += system.getSolution();
      return result;
    });
    ASSERT_TRUE(newton.getStepPolicy());
    newton.solve(state);
    using Reason = Newton::ConvergedReason;
    if (mode == 3)
    {
      EXPECT_TRUE(newton.converged());
      EXPECT_GT(calls, 1);
      EXPECT_NEAR(state[0], std::sqrt(2.0), 1e-12);
    }
    else
    {
      EXPECT_EQ(calls, 1);
      EXPECT_EQ(state[0], 1);
      EXPECT_EQ(newton.converged(), mode == 1);
      EXPECT_EQ(newton.getReport().reason,
        mode == 0     ? Reason::StepRejected
          : mode == 1 ? Reason::UserCriterion
                      : Reason::StepNormIsNotFinite);
    }
    newton.setStepPolicy(std::nullopt);
    EXPECT_FALSE(newton.getStepPolicy());
  }
}

/// @brief Clones preserve controls and solver binding, and subsequent solves reset diagnostics.
TEST(NewtonSolverTest, CloneAndReportReset)
{
  Math::Vector<Real> state(1);
  state << 1;
  ScalarNonlinearProblem problem(state);
  DenseLinearSolver linear(problem);
  Solver::NewtonSolver original(linear);
  original.setDampingFactor(0.5).setMaxIterations(1);
  original.solve(state);
  std::unique_ptr<decltype(original)> clone(original.copy());
  EXPECT_EQ(&clone->getLinearSolver(), &linear);
  const auto& constant = *clone;
  EXPECT_EQ(&constant.getLinearSolver(), &linear);
  EXPECT_EQ(clone->getDampingFactor(), 0.5);
  EXPECT_EQ(clone->getMaxIterations(), 1);
  clone->setDampingFactor(1).setMaxIterations(30).setRelativeTolerance(0);
  clone->setMonitor(std::nullopt);
  EXPECT_FALSE(clone->getMonitor());
  clone->solve(state);
  EXPECT_NEAR(state[0], std::sqrt(2.0), 1e-12);
  ASSERT_TRUE(clone->converged());
  clone->solve(state);
  EXPECT_EQ(clone->getReport().iterations, 0);
  EXPECT_EQ(clone->getReport().finalStepNorm, 0);
  EXPECT_EQ(
    clone->getReport().reason, decltype(original)::ConvergedReason::AbsoluteTolerance);
}
