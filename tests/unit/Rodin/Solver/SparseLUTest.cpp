/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::Solver
{
  using SparseLULinearSystem =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  TEST(Rodin_Solver_SparseLU, SolvesRawLinearSystem)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    SparseLULinearSystem system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 4.0;
    system.getOperator().insert(1, 0) = 2.0;
    system.getOperator().insert(0, 1) = 1.0;
    system.getOperator().insert(1, 1) = 5.0;
    system.getOperator().insert(2, 1) = 3.0;
    system.getOperator().insert(1, 2) = 1.0;
    system.getOperator().insert(2, 2) = 6.0;
    system.getVector().resize(3);
    system.getVector() << 6.0, 15.0, 24.0;

    Rodin::Solver::SparseLU solver(problem);
    EXPECT_EQ(solver.getRefinementSteps(), 0);
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getInfo().factorization, Rodin::Solver::Factorization::Numeric);

    Math::Vector<Real> expected(3);
    expected << 1.0, 2.0, 3.0;
    EXPECT_NEAR((system.getSolution() - expected).norm(), 0.0, 1e-12);
  }

  /// @brief A failed factorization is reported instead of being solved with.
  ///
  /// The solve must not run on a factorization that was never built.
  TEST(Rodin_Solver_SparseLU, ReportsFailedFactorization)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {2});
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem problem(u, v);

    // Column one is empty, so the matrix is structurally singular.
    SparseLULinearSystem system;
    system.getOperator().resize(3, 3);
    system.getOperator().insert(0, 0) = 2.0;
    system.getOperator().insert(2, 2) = 4.0;
    system.getVector().resize(3);
    system.getVector() << 1.0, 1.0, 1.0;

    // A sentinel the solver must not overwrite with a bogus solution.
    system.getSolution().resize(3);
    system.getSolution() << 7.0, 7.0, 7.0;

    Rodin::Solver::SparseLU solver(problem);
    solver.setRefinementSteps(2);
    solver.solve(system);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(system.getSolution()(0), 7.0);
    EXPECT_EQ(system.getSolution()(1), 7.0);
    EXPECT_EQ(system.getSolution()(2), 7.0);
    EXPECT_FALSE(solver.getInfo().factorization.has_value());
    EXPECT_NE(solver.getInfo().status, 0);
    EXPECT_FALSE(solver.getLastErrorMessage().empty());
  }

  template <class Scalar>
  class Rodin_Solver_SparseLURefinement : public ::testing::Test
  {
    protected:
      Rodin_Solver_SparseLURefinement()
        : m_mesh(LocalMesh::UniformGrid(Polytope::Type::Segment,{2}))
      {}
      LocalMesh m_mesh;
  };
  using RefinementScalars = ::testing::Types<Real,Complex>;
  TYPED_TEST_SUITE(Rodin_Solver_SparseLURefinement,RefinementScalars);

  /** @brief Execute exactly the configured same-precision correction updates.
   * @par Mathematical construction
   * A rank-one matrix plus a positive dyadic diagonal is nonsingular. Its
   * small diagonal amplifies factorization roundoff. The complex case has a
   * positive-definite Hermitian part. Matrix entries and the prescribed
   * integer solution give an exactly representable right-hand side.
   * An independent Eigen factorization supplies the direct solution and
   * two explicit residual-correction updates. Zero steps preserve the
   * direct result; positive counts must perform the stated updates.
   */
  TYPED_TEST(Rodin_Solver_SparseLURefinement,MatchesExplicitCorrections)
  {
    using Scalar = TypeParam;
    using System = Math::LinearSystem<Math::SparseMatrix<Scalar>,Math::Vector<Scalar>>;
    P1<Scalar> space(this->m_mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem problem(u,v);
    constexpr Eigen::Index Size = 5;
    // A small exactly representable diagonal, not an acceptance tolerance.
    const Real diagonalScale = std::ldexp(Real(1),-24);
    Scalar background(1), diagonal(1);
    if constexpr (std::is_same_v<Scalar,Complex>)
    {
      background = Scalar(1,1);
      diagonal = Scalar(1,-1);
    }
    System system;
    auto& matrix = system.getOperator();
    matrix.resize(Size,Size);
    Math::Vector<Scalar> exact(Size);
    exact << Scalar(1),Scalar(-2),Scalar(3),Scalar(-4),Scalar(5);
    for (Eigen::Index i = 0; i < Size; ++i)
      for (Eigen::Index j = 0; j < Size; ++j)
        matrix.insert(i,j) = background +
          (i == j ? diagonal * diagonalScale * Real(i+1) : Scalar(0));
    matrix.makeCompressed();
    system.getVector().resize(Size);
    for (Eigen::Index i = 0; i < Size; ++i)
      system.getVector()(i) = background*exact.sum() + diagonal*diagonalScale*Real(i+1)*exact(i);
    Eigen::SparseLU<Math::SparseMatrix<Scalar>> reference;
    reference.compute(matrix);
    ASSERT_EQ(reference.info(),Eigen::Success);
    Math::Vector<Scalar> expected = reference.solve(system.getVector());
    const Math::Vector<Scalar> direct = expected;
    Rodin::Solver::SparseLU solver(problem);
    ASSERT_EQ(solver.getRefinementSteps(),0);
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    for (Eigen::Index i = 0; i < Size; ++i)
      EXPECT_EQ(system.getSolution()(i),direct(i));
    for (size_t steps = 1; steps <= 2; ++steps)
    {
      const Math::Vector<Scalar> residual = system.getVector()-matrix*expected;
      const Math::Vector<Scalar> correction = reference.solve(residual);
      expected += correction;
      EXPECT_EQ(&solver.setRefinementSteps(steps),&solver);
      EXPECT_EQ(solver.getRefinementSteps(),steps);
      solver.solve(system);
      ASSERT_TRUE(solver.success());
      EXPECT_EQ(solver.getInfo().factorization,Rodin::Solver::Factorization::Numeric);
      for (Eigen::Index i = 0; i < Size; ++i)
        EXPECT_EQ(system.getSolution()(i),expected(i));
    }
    EXPECT_NE((expected-direct).norm(),0);
  }

  TYPED_TEST(Rodin_Solver_SparseLURefinement,CopiesConfigurationWithoutFactors)
  {
    P1<TypeParam> space(this->m_mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem problem(u,v);
    Rodin::Solver::SparseLU solver(problem);
    solver.setRefinementSteps(2);
    using System = Math::LinearSystem<Math::SparseMatrix<TypeParam>,Math::Vector<TypeParam>>;
    System system;
    system.getOperator().resize(2,2);
    system.getOperator().insert(0,0) = TypeParam(2);
    system.getOperator().insert(1,1) = TypeParam(3);
    system.getOperator().makeCompressed();
    system.getVector().resize(2);
    system.getVector() << TypeParam(6),TypeParam(12);
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    ASSERT_EQ(solver.getInfo().factorization,Rodin::Solver::Factorization::Numeric);
    auto copied(solver);
    EXPECT_EQ(copied.getRefinementSteps(),2);
    EXPECT_FALSE(copied.getInfo().factorization.has_value());
    auto moved(std::move(solver));
    EXPECT_EQ(moved.getRefinementSteps(),2);
    // No operation has failed; the existing neutral Info state is preserved.
    EXPECT_EQ(moved.success(),Rodin::Solver::Info{}.success);
    EXPECT_FALSE(moved.getInfo().factorization.has_value());
  }
}
