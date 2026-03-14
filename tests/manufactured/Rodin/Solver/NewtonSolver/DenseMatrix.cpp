/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Types.h"

using namespace Rodin;

namespace Rodin::Tests::Manufactured
{
  struct DenseLinearSolver
  {
    using MatrixType = Math::Matrix<Real>;
    using OperatorType = MatrixType;
    using VectorType = Math::Vector<Real>;
    using LinearSystemType = Math::LinearSystem<MatrixType, VectorType>;

    void solve(LinearSystemType& system)
    {
      system.getSolution() = system.getOperator().fullPivLu().solve(system.getVector());
    }
  };

  TEST(ManufacturedNewtonSolverDenseMatrix, RecoversManufacturedRoot)
  {
    Math::Vector<Real> xStar(2);
    xStar << 1.0, -1.0;

    Math::Matrix<Real> A(2, 2);
    A << 4.0, 1.0,
         1.0, 3.0;
    const Math::Vector<Real> b = A * xStar;

    Solver::NewtonSolver<Math::Vector<Real>, Math::Vector<Real>, Math::Matrix<Real>, DenseLinearSolver>
      solver(DenseLinearSolver{});

    /*
     * Residual:
     *   F_1(x) = (Ax - b)_1 + 0.2 (x_1^2 - 1)
     *   F_2(x) = (Ax - b)_2 + 0.15 (x_2^2 - 1)
     *
     * Jacobian:
     *   J(x) = A + diag(0.4 x_1, 0.3 x_2)
     *
     * Initial guess:
     *   x^(0) = [0.9, -0.8]^T
     */
    solver
      .setFunction(
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

    Math::Vector<Real> x(2);
    x << 0.9, -0.8;
    solver.solve(x);

    EXPECT_NEAR(x(0), xStar(0), 1e-10);
    EXPECT_NEAR(x(1), xStar(1), 1e-10);
  }

  TEST(ManufacturedNewtonSolverDenseMatrix, SolvesTrulyNonlinearSystem)
  {
    auto solver = Solver::NewtonSolver(DenseLinearSolver{});

    /*
     * Residual:
     *   F_1(x) = sin(x_1) + x_2 - sin(1)
     *   F_2(x) = x_1^2 + x_2^2 - 1
     *
     * Jacobian:
     *   J(x) = [ cos(x_1)   1     ]
     *          [ 2 x_1      2 x_2 ]
     *
     * Initial guess:
     *   x^(0) = [0.7, 0.3]^T
     */
    solver
      .setFunction(
        [](Math::Vector<Real>& residual, const Math::Vector<Real>& x)
        {
          residual.resize(2);
          residual(0) = std::sin(x(0)) + x(1) - std::sin(1.0);
          residual(1) = x(0) * x(0) + x(1) * x(1) - 1.0;
        })
      .setJacobian(
        [](Math::Matrix<Real>& J, const Math::Vector<Real>& x)
        {
          J.resize(2, 2);
          J(0, 0) = std::cos(x(0));
          J(0, 1) = 1.0;
          J(1, 0) = 2.0 * x(0);
          J(1, 1) = 2.0 * x(1);
        })
      .setMaxIterations(40)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-12);

    Math::Vector<Real> x(2);
    x << 0.7, 0.3;
    solver.solve(x);

    EXPECT_NEAR(x(0), 1.0, 1e-10);
    EXPECT_NEAR(x(1), 0.0, 1e-10);
  }
}
