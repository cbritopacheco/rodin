/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Types.h"

using namespace Rodin;

namespace
{
  struct DenseLinearSolver
  {
    using MatrixType = Math::Matrix<Real>;
    using VectorType = Math::Vector<Real>;
    using LinearSystemType = Math::LinearSystem<MatrixType, VectorType>;

    void solve(LinearSystemType& system)
    {
      system.getSolution() = system.getOperator().fullPivLu().solve(system.getVector());
    }
  };

  Math::Vector<Real> expectedRoot()
  {
    Math::Vector<Real> root(2);
    root << 1.0, -1.0;
    return root;
  }

  Math::Matrix<Real> baseOperator()
  {
    Math::Matrix<Real> A(2, 2);
    A << 4.0, 1.0,
         1.0, 3.0;
    return A;
  }
}

TEST(NewtonSolverTest, SolvesManufacturedDenseSystem)
{
  const auto xStar = expectedRoot();
  const auto A = baseOperator();
  const Math::Vector<Real> b = A * xStar;

  Solver::NewtonSolver<Math::Vector<Real>, Math::Vector<Real>, Math::Matrix<Real>, DenseLinearSolver>
    solver(DenseLinearSolver{});

  solver.setFunction(
      [A, b](Math::Vector<Real>& residual, const Math::Vector<Real>& x)
      {
        residual = A * x - b;
        residual(0) += 0.2 * (x(0) * x(0) - 1.0);
        residual(1) += 0.15 * (x(1) * x(1) - 1.0);
      })
    .setJacobian(
      [A](Math::Matrix<Real>& J, const Math::Vector<Real>& x)
      {
        J = A;
        J(0, 0) += 0.4 * x(0);
        J(1, 1) += 0.3 * x(1);
      })
    .setMaxIterations(30)
    .setAbsoluteTolerance(1e-12)
    .setRelativeTolerance(1e-12);

  Math::Vector<Real> x0(2);
  x0 << 0.8, -0.8;
  solver.solve(x0);

  EXPECT_NEAR(x0(0), xStar(0), 1e-10);
  EXPECT_NEAR(x0(1), xStar(1), 1e-10);
}

TEST(NewtonSolverTest, ThrowsWhenResidualNotSet)
{
  Solver::NewtonSolver<Math::Vector<Real>, Math::Vector<Real>, Math::Matrix<Real>, DenseLinearSolver>
    solver(DenseLinearSolver{});

  solver.setJacobian(
    [](Math::Matrix<Real>& J, const Math::Vector<Real>&)
    {
      J = Math::Matrix<Real>::Identity(1, 1);
    });

  Math::Vector<Real> x(1);
  x << 0.0;
  EXPECT_ANY_THROW(solver.solve(x));
}

TEST(NewtonSolverTest, ThrowsWhenJacobianNotSet)
{
  Solver::NewtonSolver<Math::Vector<Real>, Math::Vector<Real>, Math::Matrix<Real>, DenseLinearSolver>
    solver(DenseLinearSolver{});

  solver.setFunction(
    [](Math::Vector<Real>& residual, const Math::Vector<Real>& x)
    {
      residual = x;
    });

  Math::Vector<Real> x(1);
  x << 1.0;
  EXPECT_ANY_THROW(solver.solve(x));
}
