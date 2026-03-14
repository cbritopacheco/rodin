/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>

#include <gtest/gtest.h>
#include <Eigen/SparseCholesky>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Types.h"

using namespace Rodin;

namespace Rodin::Tests::Manufactured
{
  struct SparseLinearSolver
  {
    using MatrixType = Math::SparseMatrix<Real>;
    using VectorType = Math::Vector<Real>;
    using LinearSystemType = Math::LinearSystem<MatrixType, VectorType>;

    void solve(LinearSystemType& system)
    {
      Eigen::SimplicialLDLT<MatrixType> linearSolver;
      linearSolver.compute(system.getOperator());
      system.getSolution() = linearSolver.solve(system.getVector());
    }
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

    Solver::NewtonSolver<Math::Vector<Real>, Math::Vector<Real>, Math::SparseMatrix<Real>, SparseLinearSolver>
      solver(SparseLinearSolver{});

    solver
      .setFunction(
        [A, b](Math::Vector<Real>& residual, const Math::Vector<Real>& x)
        {
          residual = A * x - b;
          residual(0) += 0.2 * (x(0) * x(0) - 1.0);
          residual(1) += 0.15 * (x(1) * x(1) - 1.0);
        })
      .setJacobian(
        [A](Math::SparseMatrix<Real>& J, const Math::Vector<Real>& x)
        {
          J = A;
          J.coeffRef(0, 0) += 0.4 * x(0);
          J.coeffRef(1, 1) += 0.3 * x(1);
          J.makeCompressed();
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
}
