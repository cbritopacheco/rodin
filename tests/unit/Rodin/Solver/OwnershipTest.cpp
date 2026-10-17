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
