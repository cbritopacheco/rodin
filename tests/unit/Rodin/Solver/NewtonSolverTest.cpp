/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
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

      explicit DenseLinearSolver(ProblemBaseType& pb)
        : Parent(pb)
      {}

      void solve(DenseLinearSystem& system) override
      {
        system.getSolution() = system.getOperator().fullPivLu().solve(system.getVector());
      }

      DenseLinearSolver* copy() const noexcept override
      {
        return new DenseLinearSolver(*this);
      }
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

TEST(NewtonSolverTest, SingleTemplateParameterDeducesLinearSystem)
{
  using SolverType = Solver::NewtonSolver<DenseLinearSolver>;
  static_assert(std::is_same_v<SolverType::LinearSystemType, DenseLinearSystem>);
  SUCCEED();
}

TEST(NewtonSolverTest, PropagatesAssemblyFailure)
{
  FailingAssembleProblem pb;
  DenseLinearSolver linearSolver(pb);
  Solver::NewtonSolver newton(linearSolver);

  Math::Vector<Real> u(1);
  u << 1.0;
  EXPECT_THROW(newton.solve(u), std::runtime_error);
}
