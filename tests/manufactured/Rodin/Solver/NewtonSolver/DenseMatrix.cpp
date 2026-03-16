/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <functional>
#include <cmath>
#include <utility>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Variational/Problem.h"

using namespace Rodin;

namespace Rodin::Tests::Manufactured
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
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<Tests::Manufactured::DenseLinearSolver>
  {
    using LinearSystemType = Tests::Manufactured::DenseLinearSystem;
  };
}

namespace Rodin::Tests::Manufactured
{

  class ManufacturedDenseProblem final : public Variational::ProblemBase<DenseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      ManufacturedDenseProblem(
          Math::Vector<Real>& x,
          const Math::Matrix<Real>& A,
          const Math::Vector<Real>& b)
        : m_solution(x), m_A(A), m_b(b)
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

      ManufacturedDenseProblem& assemble() override
      {
        const auto& x = m_solution.get();
        auto& J = m_system.getOperator();
        auto& F = m_system.getVector();
        auto& du = m_system.getSolution();

        J = m_A;
        J(0, 0) += 0.4 * x(0);
        J(1, 1) += 0.3 * x(1);

        F = m_b - m_A * x;
        F(0) -= 0.2 * (x(0) * x(0) - 1.0);
        F(1) -= 0.15 * (x(1) * x(1) - 1.0);

        du.resize(2);
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

      ManufacturedDenseProblem* copy() const noexcept override
      {
        return new ManufacturedDenseProblem(*this);
      }

    private:
      std::reference_wrapper<Math::Vector<Real>> m_solution;
      Math::Matrix<Real> m_A;
      Math::Vector<Real> m_b;
      DenseLinearSystem m_system;
  };

  class StrongNonlinearDenseProblem final : public Variational::ProblemBase<DenseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      explicit StrongNonlinearDenseProblem(Math::Vector<Real>& x)
        : m_solution(x)
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

      StrongNonlinearDenseProblem& assemble() override
      {
        const auto& x = m_solution.get();
        auto& J = m_system.getOperator();
        auto& F = m_system.getVector();
        auto& du = m_system.getSolution();

        J.resize(2, 2);
        J(0, 0) = std::cos(x(0));
        J(0, 1) = 1.0;
        J(1, 0) = 2.0 * x(0);
        J(1, 1) = 2.0 * x(1);

        F.resize(2);
        F(0) = std::sin(1.0) - std::sin(x(0)) - x(1);
        F(1) = 1.0 - x(0) * x(0) - x(1) * x(1);

        du.resize(2);
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

      StrongNonlinearDenseProblem* copy() const noexcept override
      {
        return new StrongNonlinearDenseProblem(*this);
      }

    private:
      std::reference_wrapper<Math::Vector<Real>> m_solution;
      DenseLinearSystem m_system;
  };

  TEST(ManufacturedNewtonSolverDenseMatrix, RecoversManufacturedRoot)
  {
    Math::Vector<Real> xStar(2);
    xStar << 1.0, -1.0;

    Math::Matrix<Real> A(2, 2);
    A << 4.0, 1.0,
         1.0, 3.0;
    const Math::Vector<Real> b = A * xStar;

    Math::Vector<Real> x(2);
    x << 0.9, -0.8;
    ManufacturedDenseProblem pb(x, A, b);

    DenseLinearSolver linearSolver(pb);
    Solver::NewtonSolver newton(linearSolver);
    newton.setMaxIterations(30)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-12);

    newton.solve(x);

    EXPECT_NEAR(x(0), xStar(0), 1e-10);
    EXPECT_NEAR(x(1), xStar(1), 1e-10);
  }

  TEST(ManufacturedNewtonSolverDenseMatrix, SolvesStrongNonlinearSystem)
  {
    Math::Vector<Real> x(2);
    x << 0.7, 0.3;
    StrongNonlinearDenseProblem pb(x);

    DenseLinearSolver linearSolver(pb);
    Solver::NewtonSolver newton(linearSolver);
    newton.setMaxIterations(40)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-12);

    newton.solve(x);

    EXPECT_NEAR(x(0), 1.0, 1e-10);
    EXPECT_NEAR(x(1), 0.0, 1e-10);
  }
}
