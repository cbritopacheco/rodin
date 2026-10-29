/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <memory>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/GMRES.h"
#include "Rodin/Variational.h"
#include "../../../benchmarks/Workloads/StokesBlockPreconditioner.h"

namespace Rodin::Tests::Unit
{
  /// @brief Configurable test preconditioner with observable setup and failure.
  template <class Scalar>
  class CountingDiagonal
  {
    public:
      using StorageIndex = int;
      enum
      {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic
      };
      Real scale = 1;
      bool fail = false;
      std::shared_ptr<size_t> setups = std::make_shared<size_t>(0);

      CountingDiagonal() = default;
      CountingDiagonal& operator=(const CountingDiagonal& other)
      {
        scale = other.scale;
        fail = other.fail;
        setups = other.setups;
        m_status = Eigen::InvalidInput;
        return *this;
      }

      template <class Matrix>
      CountingDiagonal& compute(const Matrix& matrix)
      {
        ++*setups;
        m_status = fail ? Eigen::NumericalIssue : Eigen::Success;
        if (!fail)
          m_diagonal.compute(matrix);
        return *this;
      }

      template <class Rhs, class Destination>
      void _solve_impl(const Rhs& rhs, Destination& result) const
      {
        result = scale * m_diagonal.solve(rhs);
      }

      template <class Rhs>
      const Eigen::Solve<CountingDiagonal, Rhs> solve(
        const Eigen::MatrixBase<Rhs>& rhs) const
      {
        return Eigen::Solve<CountingDiagonal, Rhs>(*this, rhs.derived());
      }

      Eigen::ComputationInfo info() const
      {
        return m_status;
      }
      Eigen::Index rows() const
      {
        return m_diagonal.rows();
      }
      Eigen::Index cols() const
      {
        return m_diagonal.cols();
      }

    private:
      Eigen::DiagonalPreconditioner<Scalar> m_diagonal;
      Eigen::ComputationInfo m_status = Eigen::InvalidInput;
  };

  /// @brief Verifies configuration ownership, backend setup failure and recovery.
  template <class Scalar, class Matrix>
  void checkCustomGMRES()
  {
    using namespace Variational;
    using System = Math::LinearSystem<Matrix, Math::Vector<Scalar>>;
    using SolverType = Solver::GMRES<System, CountingDiagonal<Scalar>>;
    auto mesh = Geometry::LocalMesh::UniformGrid(Geometry::Polytope::Type::Segment, {8});
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
    for (Eigen::Index i = 0; i < matrix.rows(); ++i)
    {
      matrix.coeffRef(i, i) = Scalar(4);
      exact(i) = Scalar(i + 1);
      if (i + 1 < matrix.cols())
        matrix.coeffRef(i, i + 1) = Scalar(0.5);
    }
    system.getVector() = matrix * exact;
    SolverType solver(problem);
    solver.setTolerance(1e-13).setRestart(4).setMaxIterations(200);
    solver.getPreconditioner().scale = 0.5;
    const auto setups = solver.getPreconditioner().setups;
    const auto check = [&](auto& tested) {
      const auto before = *setups;
      system.getSolution().setZero();
      tested.solve(system);
      EXPECT_TRUE(tested.success());
      EXPECT_GT(*setups, before);
      EXPECT_LT((system.getSolution() - exact).norm() / exact.norm(), 1e-11);
      EXPECT_LT((matrix * system.getSolution() - system.getVector()).norm(), 1e-10);
      EXPECT_EQ(std::as_const(tested).getPreconditioner().scale, 0.5);
    };
    check(solver);
    SolverType copied(solver);
    check(copied);
    SolverType moved(std::move(copied));
    check(moved);
    std::unique_ptr<SolverType> cloned(solver.copy());
    check(*cloned);
    solver.getPreconditioner().fail = true;
    system.getSolution().setConstant(Scalar(-3));
    const auto sentinel = system.getSolution().eval();
    solver.solve(system);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ((system.getSolution() - sentinel).norm(), 0);
    solver.getPreconditioner().fail = false;
    check(solver);
  }

  TEST(Rodin_Solver_GMRESPreconditioner, SparseReal)
  {
    checkCustomGMRES<Real, Math::SparseMatrix<Real>>();
  }
  TEST(Rodin_Solver_GMRESPreconditioner, SparseComplex)
  {
    checkCustomGMRES<Complex, Math::SparseMatrix<Complex>>();
  }
  TEST(Rodin_Solver_GMRESPreconditioner, DenseReal)
  {
    checkCustomGMRES<Real, Math::Matrix<Real>>();
  }
  TEST(Rodin_Solver_GMRESPreconditioner, DenseComplex)
  {
    checkCustomGMRES<Complex, Math::Matrix<Complex>>();
  }

  /// @brief Builds a small mixed system whose Schur complement is known exactly.
  Math::Matrix<Real> blockOperator()
  {
    Math::Matrix<Real> matrix = Math::Matrix<Real>::Zero(6, 6);
    matrix(0, 0) = 2;
    matrix(1, 1) = 3;
    matrix(2, 2) = 4;
    matrix(3, 0) = 1;
    matrix(3, 1) = 2;
    matrix(4, 1) = 1;
    matrix(4, 2) = 1;
    matrix.block(0, 3, 3, 2) = -matrix.block(3, 0, 2, 3).transpose().eval();
    matrix(3, 5) = matrix(5, 3) = 1;
    matrix(4, 5) = matrix(5, 4) = 2;
    return matrix;
  }

  /// @brief Returns the exact pressure Schur complement of blockOperator.
  Math::SparseMatrix<Real> blockMass()
  {
    const auto matrix = blockOperator();
    const Math::Matrix<Real> mass = matrix.block(3, 0, 2, 3) *
      matrix.topLeftCorner(3, 3).diagonal().cwiseInverse().asDiagonal() *
      matrix.block(3, 0, 2, 3).transpose();
    return mass.sparseView();
  }

  TEST(Rodin_Solver_StokesBlockPreconditioner, BlockSignsAndNonzeroMeanConstraint)
  {
    const auto dense = blockOperator();
    const Math::SparseMatrix<Real> matrix = dense.sparseView();
    Benchmarks::StokesBlockPreconditioner preconditioner;
    preconditioner.configure(blockMass(), 3);
    preconditioner.compute(matrix);
    ASSERT_EQ(preconditioner.info(), Eigen::Success);
    for (Eigen::Index j = 0; j < dense.cols(); ++j)
    {
      Math::Vector<Real> exact = Math::Vector<Real>::Zero(6);
      exact(j) = 1;
      const Math::Vector<Real> rhs = dense * exact;
      const Math::Vector<Real> result = preconditioner.solve(rhs);
      EXPECT_LT((result - exact).norm(), 1e-11);
      EXPECT_NEAR(result(3) + 2 * result(4), rhs(5), 1e-12);
    }
  }

  TEST(
    Rodin_Solver_StokesBlockPreconditioner, ConfigurationOutlivesSourceAndRebuildsFactors)
  {
    std::unique_ptr<Benchmarks::StokesBlockPreconditioner> copied;
    const Math::SparseMatrix<Real> matrix = blockOperator().sparseView();
    {
      Benchmarks::StokesBlockPreconditioner original;
      original.configure(blockMass(), 3);
      original.compute(matrix);
      ASSERT_EQ(original.info(), Eigen::Success);
      copied = std::make_unique<Benchmarks::StokesBlockPreconditioner>(original);
    }
    EXPECT_EQ(copied->info(), Eigen::InvalidInput);
    copied->compute(matrix);
    ASSERT_EQ(copied->info(), Eigen::Success);
    const Math::Vector<Real> rhs = Math::Vector<Real>::Ones(6);
    const Math::Vector<Real> solution = copied->solve(rhs);
    EXPECT_LT((matrix * solution - rhs).norm(), 1e-11);
  }

  TEST(
    Rodin_Solver_StokesBlockPreconditioner, InvalidConfigurationAndSingularMassFailSetup)
  {
    const Math::SparseMatrix<Real> matrix = blockOperator().sparseView();
    Benchmarks::StokesBlockPreconditioner preconditioner;
    EXPECT_EQ(preconditioner.compute(matrix).info(), Eigen::InvalidInput);
    preconditioner.configure(blockMass(), 2);
    EXPECT_EQ(preconditioner.compute(matrix).info(), Eigen::InvalidInput);
    preconditioner.configure(Math::SparseMatrix<Real>(2, 2), 3);
    EXPECT_EQ(preconditioner.compute(matrix).info(), Eigen::NumericalIssue);
  }
}
