/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Solver ownership and refactorization against independent algebraic systems. */
#include <gtest/gtest.h>
#include <memory>
#include <utility>
#include <type_traits>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational.h"
#include RODIN_TEST_SOLVER_HEADER

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar, bool Sparse>
  void checkOwnership()
  {
    using Operator =
      std::conditional_t<Sparse, Math::SparseMatrix<Scalar>, Math::Matrix<Scalar>>;
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    using Algorithm = Rodin::Solver::RODIN_TEST_SOLVER<System>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem<System, decltype(u), decltype(v)> problem(u, v);
    std::unique_ptr<Algorithm> survivor;
    System system;
    Math::Matrix<Scalar> matrix(4, 4);
    Math::Vector<Scalar> expected(4);
    constexpr Real SolutionTolerance = 1e-10;
    constexpr Real ResidualTolerance = 1e-11;
    // A full four-dimensional Krylov cycle needs four steps in exact arithmetic.
    // Twenty-five such cycles are a heuristic work budget; accuracy is checked separately.
    constexpr size_t RestartDimension = 4;
    constexpr size_t IterationLimit = 25 * RestartDimension;
    // Exercise DGMRES recycling with one retained direction on this small system.
    constexpr size_t DeflationVectorCount = 1;
    auto check = [&](auto& candidate, size_t phase) {
      SCOPED_TRACE(phase);
      EXPECT_EQ(&candidate.getProblem(), &problem);
      const auto& readOnly = candidate;
      EXPECT_EQ(&readOnly.getProblem(), &problem);
      // Strict diagonal dominance gives positive definiteness in every phase.
      matrix.setZero();
      for (Eigen::Index i = 0; i < 4; ++i)
      {
        matrix(i, i) = Scalar(6 + i + phase);
        expected[i] = Scalar((i % 2 ? -1 : 1) * Real(i + phase + 1));
      }
      Scalar offDiagonal = Scalar(1);
      if constexpr (std::is_same_v<Scalar, Complex>)
      {
        offDiagonal += Scalar(0, 0.5);
        expected[1] += Scalar(0, Real(phase + 1));
      }
      for (Eigen::Index i = 0; i < 3; ++i)
      {
        matrix(i, i + 1) = offDiagonal;
        matrix(i + 1, i) = Eigen::numext::conj(offDiagonal);
      }
      if constexpr (Sparse)
        system.getOperator() = matrix.sparseView();
      else
        system.getOperator() = matrix;
      system.getVector() = matrix * expected;
      // Exercise both a cold solve and a correctly sized, deliberately wrong guess.
      if (phase % 2)
        system.getSolution().setConstant(4, Scalar(7));
      else
        system.getSolution().resize(0);
      if constexpr (requires { candidate.setTolerance(ResidualTolerance); })
        candidate.setTolerance(ResidualTolerance).setMaxIterations(IterationLimit);
      if constexpr (requires { candidate.setRestart(RestartDimension); })
        candidate.setRestart(RestartDimension);
      if constexpr (requires { candidate.setEigenv(DeflationVectorCount); })
        candidate.setEigenv(DeflationVectorCount);
      candidate.solve(system);
      if constexpr (requires { candidate.success(); })
      {
        EXPECT_TRUE(candidate.success());
      }
      ASSERT_EQ(system.getSolution().size(), expected.size());
      EXPECT_LT((system.getSolution() - expected).norm(), SolutionTolerance);
      EXPECT_LT((matrix * system.getSolution() - system.getVector()).norm() /
          system.getVector().norm(),
        ResidualTolerance);
    };
    {
      Algorithm original(problem);
      check(original, 0);
      Algorithm copied(original);
      check(copied, 1);
      Algorithm moved(std::move(copied));
      check(moved, 2);
      const Rodin::Solver::LinearSolverBase<System>& polymorphic = original;
      survivor.reset(static_cast<Algorithm*>(polymorphic.copy()));
      check(*survivor, 3);
      check(original, 4);
    }
    // The clone retains the borrowed problem, not the original solver's factors.
    check(*survivor, 5);
  }

#ifdef RODIN_TEST_KRYLOV_CONFIGURATION
  /**
   * @brief Checks configured cold and previously used solvers against Eigen after ownership changes.
   *
   * Settings are applied only to the original. The one-step budget exposes lost
   * iteration limits, the loose tolerance exposes lost tolerances, and restarted
   * solves exercise the requested Krylov and deflation settings. Default limits
   * remain size-dependent when the source has not been given a matrix.
   *
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Sparse Whether the operator uses sparse storage.
   */
  template <class Scalar, bool Sparse>
  void checkConfiguration()
  {
    using Operator =
      std::conditional_t<Sparse, Math::SparseMatrix<Scalar>, Math::Matrix<Scalar>>;
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    using Algorithm = Rodin::Solver::RODIN_TEST_SOLVER<System>;
    using Reference = Eigen::RODIN_TEST_SOLVER<Operator>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem<System, decltype(u), decltype(v)> problem(u, v);
    constexpr Eigen::Index Size = 8;
    Math::Matrix<Scalar> matrix = Math::Matrix<Scalar>::Zero(Size, Size);
    Math::Vector<Scalar> exact(Size);
    for (Eigen::Index i = 0; i < Size; ++i)
    {
      matrix(i, i) = Scalar(4 + i);
      exact[i] = Scalar(i + 1);
      if constexpr (std::is_same_v<Scalar, Complex>)
        exact[i] += Scalar(0, Real(i % 3 + 1) / 2);
      if (i + 1 < Size)
      {
        matrix(i, i + 1) = Scalar(1);
        matrix(i + 1, i) = Scalar(-0.5);
      }
    }
    for (const size_t configuration : {size_t(0), size_t(1), size_t(2), size_t(3)})
    {
      SCOPED_TRACE(configuration);
      const Real tolerance = configuration == 1 ? Real(0.25) : Real(1e-13);
      const size_t budget = configuration == 0 ? 1 : 200;
      constexpr size_t Restart = 4;
      // Memcheck detects uninitialized Schur workspace in Eigen 3.4 DGMRES
      // for both fields. Nonzero settings are copied in single-cycle solves;
      // restarted numerical comparisons disable that separate backend defect.
      const size_t deflation = configuration < 2 ? 1 : 0;
      System system;
      if constexpr (Sparse)
        system.getOperator() = matrix.sparseView();
      else
        system.getOperator() = matrix;
      system.getVector() = matrix * exact;
      Reference reference;
      reference.setTolerance(tolerance);
      if (configuration != 3)
        reference.setMaxIterations(budget);
      reference.set_restart(Restart);
      if constexpr (requires { reference.setEigenv(deflation); })
        reference.setEigenv(deflation);
      reference.compute(system.getOperator());
      const Math::Vector<Scalar> zero = Math::Vector<Scalar>::Zero(Size);
      const Math::Vector<Scalar> expected =
        reference.solveWithGuess(system.getVector(), zero);
      const bool expectedSuccess = reference.info() == Eigen::Success;
      const auto check = [&](Algorithm& solver, const char* ownership) {
        SCOPED_TRACE(ownership);
        EXPECT_EQ(&solver.getProblem(), &problem);
        system.getSolution() = zero;
        solver.solve(system);
        EXPECT_EQ(solver.success(), expectedSuccess);
        EXPECT_LT((system.getSolution() - expected).norm(), 1e-11);
        if constexpr (requires { solver.getIterationNumber(); })
        {
          EXPECT_EQ(solver.getIterationNumber(), size_t(reference.iterations()));
        }
      };
      std::unique_ptr<Algorithm> clone;
      {
        Algorithm original(problem);
        original.setTolerance(tolerance).setRestart(Restart);
        if (configuration != 3)
          original.setMaxIterations(budget);
        if constexpr (requires { original.setEigenv(deflation); })
          original.setEigenv(deflation);
        else if constexpr (requires { original.setDeflationSize(deflation); })
          original.setDeflationSize(deflation);
        // Copy the default-limit case before any matrix has been supplied.
        Algorithm coldCopy(original);
        check(coldCopy, "cold copy");
        check(original, "original");
        Algorithm copied(original);
        check(copied, "warm copy");
        Algorithm moved(std::move(copied));
        check(moved, "move");
        clone.reset(original.copy());
      }
      check(*clone, "clone");
    }
  }

  /// @brief Real dense and sparse Krylov solvers retain settings without reconfiguration.
  TEST(SolverOwnership, ConfiguredReal)
  {
    checkConfiguration<Real, false>();
    checkConfiguration<Real, true>();
  }

  /// @brief Complex dense and sparse Krylov solvers retain settings without reconfiguration.
  TEST(SolverOwnership, ConfiguredComplex)
  {
    checkConfiguration<Complex, false>();
    checkConfiguration<Complex, true>();
  }
#endif

  /// @brief Real solvers preserve their problem and recompute factors after ownership changes.
  TEST(SolverOwnership, Real)
  {
#ifdef RODIN_TEST_DENSE
    checkOwnership<Real, false>();
#endif
#ifdef RODIN_TEST_SPARSE
    checkOwnership<Real, true>();
#endif
  }

  /// @brief Complex Hermitian systems check conjugation and ownership on each supported storage.
  TEST(SolverOwnership, Complex)
  {
#ifdef RODIN_TEST_DENSE
    checkOwnership<Complex, false>();
#endif
#ifdef RODIN_TEST_SPARSE
    checkOwnership<Complex, true>();
#endif
  }
}
