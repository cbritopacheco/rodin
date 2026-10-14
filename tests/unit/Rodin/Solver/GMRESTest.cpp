/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Sparse and dense GMRES diagnostics agree with the backend's last solve.
 */
#include <gtest/gtest.h>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/GMRES.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Checks diagnostics after convergence, failure and an exact initial guess.
   *
   * A fixed-layout problem supplies the system storage. Its entries are replaced
   * by a strictly diagonally dominant nonsymmetric tridiagonal operator and a
   * load computed from a known coefficient vector. Independent Eigen solves
   * verify the reported count and estimate; true residuals verify the solution.
   * Dyadic entries and coefficients make the exact-guess residual bitwise
   * zero. This is required by Eigen's zero-initial-residual branch; a residual
   * at roundoff scale is not equivalent to exact zero for that branch.
   *
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Matrix Sparse or dense native operator representation.
   */
  template <class Scalar, class Matrix>
  void checkGMRESDiagnostics()
  {
    using namespace Variational;
    using System = Math::LinearSystem<Matrix, Math::Vector<Scalar>>;
    constexpr Real tolerance = 1e-13;
    constexpr size_t restart = 3, maximumIterations = 200;
    auto mesh =
      Geometry::Mesh<Context::Local>::UniformGrid(Geometry::Polytope::Type::Segment, {8});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem<System, decltype(u), decltype(v)> problem(u, v);
    problem = Integral(u, v);
    problem.assemble();
    auto& system = problem.getLinearSystem();
    auto& matrix = system.getOperator();
    matrix.setZero();
    Math::Vector<Scalar> exact(matrix.rows());
    const size_t rows = static_cast<size_t>(matrix.rows());
    for (size_t j = 0; j < rows; ++j)
    {
      if constexpr (std::is_same_v<Scalar, Complex>)
      {
        exact(j) = Scalar(1 + Real(j), Real(j % 3 + 1) / 2);
        matrix.coeffRef(j, j) = Scalar(4, 0.5);
      }
      else
      {
        exact(j) = 1 + Real(j);
        matrix.coeffRef(j, j) = 4;
      }
      if (j > 0)
        matrix.coeffRef(j, j - 1) = -0.5;
      if (j + 1 < rows)
        matrix.coeffRef(j, j + 1) = 1;
    }
    system.getVector() = matrix * exact;
    Solver::GMRES solver(problem);
    Eigen::GMRES<Matrix> reference;
    solver.setTolerance(tolerance).setRestart(restart);
    reference.setTolerance(tolerance);
    reference.set_restart(restart);
    for (const size_t cap : {maximumIterations, size_t(1), maximumIterations})
    {
      solver.setMaxIterations(cap);
      reference.setMaxIterations(cap);
      system.getSolution().setZero();
      reference.compute(matrix);
      const Math::Vector<Scalar> expected =
        reference.solveWithGuess(system.getVector(), system.getSolution());
      solver.solve(system);
      EXPECT_EQ(solver.getIterationNumber(), size_t(reference.iterations()));
      EXPECT_EQ(solver.getError(), reference.error());
      EXPECT_EQ(solver.success(), reference.info() == Eigen::Success);
      EXPECT_EQ((system.getSolution() - expected).norm(), 0);
      EXPECT_EQ(solver.success(), cap != 1);
      if (solver.success())
      {
        EXPECT_LT((matrix * system.getSolution() - system.getVector()).norm() /
            system.getVector().norm(),
          1e-11);
        EXPECT_LT((system.getSolution() - exact).norm() / exact.norm(), 1e-11);
      }
    }
    system.getSolution() = exact;
    const Math::Vector<Scalar> initialResidual = system.getVector() - matrix * exact;
    EXPECT_EQ(initialResidual.norm(), 0);
    solver.setMaxIterations(1).solve(system);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0);
    EXPECT_LT(solver.getError(), tolerance);
  }

  TEST(Rodin_Solver_GMRES, SparseRealDiagnostics)
  {
    checkGMRESDiagnostics<Real, Math::SparseMatrix<Real>>();
  }

  TEST(Rodin_Solver_GMRES, SparseComplexDiagnostics)
  {
    checkGMRESDiagnostics<Complex, Math::SparseMatrix<Complex>>();
  }

  TEST(Rodin_Solver_GMRES, DenseRealDiagnostics)
  {
    checkGMRESDiagnostics<Real, Math::Matrix<Real>>();
  }

  TEST(Rodin_Solver_GMRES, DenseComplexDiagnostics)
  {
    checkGMRESDiagnostics<Complex, Math::Matrix<Complex>>();
  }
}
