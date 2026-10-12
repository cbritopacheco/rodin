/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Dense/sparse iterative solver specialization contracts. */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/IDRS.h"
#include "Rodin/Solver/LeastSquaresCG.h"
#include "Rodin/Solver/MINRES.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <template <class> class Algorithm, class Scalar, bool Sparse>
  void checkIterativeSolver()
  {
    using Operator =
      std::conditional_t<Sparse, Math::SparseMatrix<Scalar>, Math::Matrix<Scalar>>;
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem<System, decltype(u), decltype(v)> problem(u, v);
    Algorithm<System> solver(problem);
    solver.setTolerance(1e-12).setMaxIterations(200);
    if constexpr (requires { solver.setS(2); })
      solver.setS(2);
    if constexpr (requires { solver.setL(2); })
      solver.setL(2);
    Math::Matrix<Scalar> matrix(4, 4);
    matrix << 5, 1, 0, 0, 1, 6, 1, 0, 0, 1, 7, 1, 0, 0, 1, 8;
    if constexpr (std::is_same_v<Scalar, Complex>)
    {
      matrix(0, 1) += Complex(0, 0.5);
      matrix(1, 0) -= Complex(0, 0.5);
    }
    Math::Vector<Scalar> expected(4);
    expected << 1, -2, 3, -4;
    if constexpr (std::is_same_v<Scalar, Complex>)
      expected(2) += Complex(0, 2);
    System system;
    if constexpr (Sparse)
      system.getOperator() = matrix.sparseView();
    else
      system.getOperator() = matrix;
    system.getVector() = matrix * expected;
    solver.solve(system);
    if constexpr (requires { solver.success(); })
      EXPECT_TRUE(solver.success());
    EXPECT_LT((system.getSolution() - expected).norm(), 1e-9);
    std::unique_ptr<Algorithm<System>> clone(solver.copy());
    system.getSolution().resize(0);
    clone->solve(system);
    EXPECT_LT((system.getSolution() - expected).norm(), 1e-9);
    Algorithm<System> moved(std::move(*clone));
    system.getSolution().resize(0);
    moved.solve(system);
    EXPECT_LT((system.getSolution() - expected).norm(), 1e-9);
    // A nonzero guess takes the separate solveWithGuess path.
    system.getSolution().setConstant(Scalar(0.25));
    solver.solve(system);
    if constexpr (requires { solver.success(); })
      EXPECT_TRUE(solver.success());
    EXPECT_LT((system.getSolution() - expected).norm(), 1e-9);
  }

  template <template <class> class Algorithm>
  void checkStorageAndScalarVariants()
  {
    checkIterativeSolver<Algorithm, Real, false>();
    checkIterativeSolver<Algorithm, Real, true>();
    checkIterativeSolver<Algorithm, Complex, false>();
    checkIterativeSolver<Algorithm, Complex, true>();
  }
}

TEST(Rodin_Solver_Iterative, MINRESDenseSparseReal)
{
  // Eigen 3.4 MINRES only implements real symmetric arithmetic.
  checkIterativeSolver<Rodin::Solver::MINRES, Real, false>();
  checkIterativeSolver<Rodin::Solver::MINRES, Real, true>();
}

TEST(Rodin_Solver_Iterative, LeastSquaresCGDenseSparseRealComplex)
{
  checkStorageAndScalarVariants<Rodin::Solver::LeastSquaresCG>();
}

TEST(Rodin_Solver_Iterative, IDRSDenseSparseRealComplex)
{
  checkStorageAndScalarVariants<Rodin::Solver::IDRS>();
}
