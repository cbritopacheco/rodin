/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief ParU input rejection, retained factors and singular-system recovery. */
#include <gtest/gtest.h>
#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/ParU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  class ParUContract : public ::testing::Test
  {
    protected:
      ParUContract()
        : mesh(LocalMesh::UniformGrid(Polytope::Type::Segment, {2})),
          space(mesh),
          u(space),
          v(space),
          problem(u, v)
      {}

      LocalMesh mesh;
      P1<Real> space;
      TrialFunction<GridFunction<P1<Real>, Math::Vector<Real>>, P1<Real>> u;
      TestFunction<P1<Real>> v;
      Problem<System, decltype(u), decltype(v)> problem;
  };

  /// @brief Nonsquare input is rejected before creating backend resources or altering data.
  TEST_F(ParUContract, NonsquareInputPreservesSystem)
  {
    System system;
    system.getOperator().resize(2, 3);
    system.getOperator().insert(0, 0) = 2;
    system.getOperator().insert(1, 2) = 4;
    system.getVector() = Math::Vector<Real>::Constant(2, 3);
    system.getSolution() = Math::Vector<Real>::Constant(2, 7);
    const Math::Matrix<Real> before(system.getOperator());
    Rodin::Solver::ParU solver(problem);
    EXPECT_THROW(solver.factorize(system), Alert::Exception);
    EXPECT_EQ((Math::Matrix<Real>(system.getOperator()) - before).norm(), 0);
    EXPECT_EQ(system.getVector().sum(), 6);
    EXPECT_EQ(system.getSolution().sum(), 14);
    EXPECT_EQ(solver.getResources().control, nullptr);
    EXPECT_EQ(solver.getResources().symbolic, nullptr);
    EXPECT_EQ(solver.getResources().numeric, nullptr);
  }

  /// @brief A wrong load size preserves numeric factors and permits a corrected solve.
  TEST_F(ParUContract, BadLoadPreservesFactorsAndSolution)
  {
    System system;
    system.getOperator().resize(2, 2);
    system.getOperator().insert(0, 0) = 2;
    system.getOperator().insert(1, 1) = 4;
    system.getVector() = Math::Vector<Real>::Constant(3, 5);
    system.getSolution() = Math::Vector<Real>::Constant(2, 7);
    Rodin::Solver::ParU solver(problem);
    solver.setMaxThreads(1).factorize(system);
    ASSERT_TRUE(solver.success());
    const auto symbolic = solver.getResources().symbolic;
    const auto numeric = solver.getResources().numeric;
    EXPECT_THROW(solver.solve(system), Alert::Exception);
    EXPECT_EQ(solver.getResources().symbolic, symbolic);
    EXPECT_EQ(solver.getResources().numeric, numeric);
    EXPECT_EQ(system.getVector().size(), 3);
    EXPECT_EQ(system.getSolution().sum(), 14);
    system.getVector().resize(2);
    system.getVector() << 2, 8;
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(system.getSolution()[0], 1);
    EXPECT_EQ(system.getSolution()[1], 2);
  }

  /// @brief Numerical singularity retains symbolic analysis and leaves the output untouched.
  TEST_F(ParUContract, SingularFactorizationRecoversAfterValueChange)
  {
    System system;
    system.getOperator().resize(2, 2);
    system.getOperator().insert(0, 0) = 1;
    system.getOperator().insert(1, 0) = 2;
    system.getOperator().insert(0, 1) = 2;
    system.getOperator().insert(1, 1) = 4;
    system.getVector() = Math::Vector<Real>::Constant(2, 1);
    system.getSolution() = Math::Vector<Real>::Constant(2, 7);
    Rodin::Solver::ParU solver(problem);
    solver.setMaxThreads(1).solve(system);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getInfo().factorization, Rodin::Solver::Factorization::Symbolic);
    EXPECT_NE(solver.getResources().symbolic, nullptr);
    EXPECT_EQ(solver.getResources().numeric, nullptr);
    EXPECT_EQ(system.getSolution().sum(), 14);
    const auto symbolic = solver.getResources().symbolic;
    system.getOperator().coeffRef(1, 1) = 5;
    system.getVector() << 0, -1;
    solver.factorize(system);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getResources().symbolic, symbolic);
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    constexpr Real SolutionTolerance = 1e-12;
    EXPECT_NEAR(system.getSolution()[0], 2, SolutionTolerance);
    EXPECT_NEAR(system.getSolution()[1], -1, SolutionTolerance);
  }
}
