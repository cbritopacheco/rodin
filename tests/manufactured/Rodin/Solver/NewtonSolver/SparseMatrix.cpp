/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>
#include <cmath>
#include <functional>
#include <utility>

#include <gtest/gtest.h>
#include <Eigen/SparseCholesky>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Variational/Problem.h"

using namespace Rodin;

namespace Rodin::Tests::Manufactured
{
  using SparseLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  class SparseLinearSolver final : public Solver::LinearSolverBase<SparseLinearSystem>
  {
    public:
      using Parent = Solver::LinearSolverBase<SparseLinearSystem>;
      using ProblemBaseType = typename Parent::ProblemBaseType;

      using Parent::solve;

      explicit SparseLinearSolver(ProblemBaseType& pb)
        : Parent(pb)
      {}

      void solve(SparseLinearSystem& system) override
      {
        Eigen::SimplicialLDLT<Math::SparseMatrix<Real>> linearSolver;
        linearSolver.compute(system.getOperator());
        system.getSolution() = linearSolver.solve(system.getVector());
      }

      SparseLinearSolver* copy() const noexcept override
      {
        return new SparseLinearSolver(*this);
      }
  };
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<Tests::Manufactured::SparseLinearSolver>
  {
    using LinearSystemType = Tests::Manufactured::SparseLinearSystem;
  };
}

namespace Rodin::Tests::Manufactured
{

  class ManufacturedSparseProblem final : public Variational::ProblemBase<SparseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<SparseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      ManufacturedSparseProblem(
          Math::Vector<Real>& x,
          const Math::SparseMatrix<Real>& A,
          const Math::Vector<Real>& b)
        : m_solution(x), m_A(A), m_b(b)
      {}

      Parent& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Solver::LinearSolverBase<SparseLinearSystem>& solver) override
      {
        assemble();
        solver.solve(m_system);
      }

      ManufacturedSparseProblem& assemble() override
      {
        const auto& x = m_solution.get();
        auto& J = m_system.getOperator();
        auto& F = m_system.getVector();
        auto& du = m_system.getSolution();

        J = m_A;
        J.coeffRef(0, 0) += 0.4 * x(0);
        J.coeffRef(1, 1) += 0.3 * x(1);
        J.makeCompressed();

        F = m_b - m_A * x;
        F(0) -= 0.2 * (x(0) * x(0) - 1.0);
        F(1) -= 0.15 * (x(1) * x(1) - 1.0);

        du.resize(2);
        return *this;
      }

      SparseLinearSystem& getLinearSystem() override
      {
        return m_system;
      }

      const SparseLinearSystem& getLinearSystem() const override
      {
        return m_system;
      }

      ManufacturedSparseProblem* copy() const noexcept override
      {
        return new ManufacturedSparseProblem(*this);
      }

    private:
      std::reference_wrapper<Math::Vector<Real>> m_solution;
      Math::SparseMatrix<Real> m_A;
      Math::Vector<Real> m_b;
      SparseLinearSystem m_system;
  };

  class StrongNonlinearSparseProblem final : public Variational::ProblemBase<SparseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<SparseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      explicit StrongNonlinearSparseProblem(Math::Vector<Real>& x)
        : m_solution(x)
      {}

      Parent& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Solver::LinearSolverBase<SparseLinearSystem>& solver) override
      {
        assemble();
        solver.solve(m_system);
      }

      StrongNonlinearSparseProblem& assemble() override
      {
        const auto& x = m_solution.get();
        auto& J = m_system.getOperator();
        auto& F = m_system.getVector();
        auto& du = m_system.getSolution();

        std::vector<Eigen::Triplet<Real>> entries;
        entries.emplace_back(0, 0, std::cos(x(0)));
        entries.emplace_back(0, 1, 1.0);
        entries.emplace_back(1, 0, 2.0 * x(0));
        entries.emplace_back(1, 1, 2.0 * x(1));
        J.resize(2, 2);
        J.setFromTriplets(entries.begin(), entries.end());

        F.resize(2);
        F(0) = std::sin(1.0) - std::sin(x(0)) - x(1);
        F(1) = 1.0 - x(0) * x(0) - x(1) * x(1);

        du.resize(2);
        return *this;
      }

      SparseLinearSystem& getLinearSystem() override
      {
        return m_system;
      }

      const SparseLinearSystem& getLinearSystem() const override
      {
        return m_system;
      }

      StrongNonlinearSparseProblem* copy() const noexcept override
      {
        return new StrongNonlinearSparseProblem(*this);
      }

    private:
      std::reference_wrapper<Math::Vector<Real>> m_solution;
      SparseLinearSystem m_system;
  };

  TEST(ManufacturedNewtonSolverSparseMatrix, RecoversManufacturedRoot)
  {
    Math::Vector<Real> xStar(2);
    xStar << 1.0, -1.0;

    Math::SparseMatrix<Real> A(2, 2);
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.emplace_back(0, 0, 4.0);
    triplets.emplace_back(0, 1, 1.0);
    triplets.emplace_back(1, 0, 1.0);
    triplets.emplace_back(1, 1, 3.0);
    A.setFromTriplets(triplets.begin(), triplets.end());

    const Math::Vector<Real> b = A * xStar;

    Math::Vector<Real> x(2);
    x << 0.9, -0.8;
    ManufacturedSparseProblem pb(x, A, b);

    SparseLinearSolver linearSolver(pb);
    Solver::NewtonSolver newton(linearSolver);
    newton.setMaxIterations(30)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-12);

    newton.solve(x);

    EXPECT_NEAR(x(0), xStar(0), 1e-10);
    EXPECT_NEAR(x(1), xStar(1), 1e-10);
  }

  TEST(ManufacturedNewtonSolverSparseMatrix, SolvesStrongNonlinearSystem)
  {
    Math::Vector<Real> x(2);
    x << 0.7, 0.3;
    StrongNonlinearSparseProblem pb(x);

    SparseLinearSolver linearSolver(pb);
    Solver::NewtonSolver newton(linearSolver);
    newton.setMaxIterations(40)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-12);

    newton.solve(x);

    EXPECT_NEAR(x(0), 1.0, 1e-10);
    EXPECT_NEAR(x(1), 0.0, 1e-10);
  }
}
