/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/// @file @brief SNES iterate identity and state-synchronization contract.

#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/PETSc.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_PETSc_SNESState, DistinctIteratesWithEqualCountersAreNotCachedTogether)
{
  auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {3});
  P1<PetscScalar> space(mesh);
  PETSc::Variational::GridFunction state(space);
  state = Zero();
  PETSc::Variational::TrialFunction du(space);
  PETSc::Variational::TestFunction v(space);
  Problem problem(du, v);
#ifdef PETSC_USE_COMPLEX
  const auto one = ComplexFunction(1);
#else
  const auto one = RealFunction(1);
#endif
  problem = Integral((one + 3 * state * state) * du, v) +
    Integral(state + state * state * state, v);
  problem.assemble();
  Solver::KSP ksp(problem);
  Solver::SNES snes(ksp);
  size_t updates = 0;
  snes.setStateUpdate([&](const PETSc::Math::Vector& x) {
    ++updates;
    state.setData(x);
  });

  ::Vec first = nullptr, second = nullptr, residual = nullptr, reference = nullptr;
  ASSERT_EQ(VecDuplicate(problem.getLinearSystem().getSolution(), &first), PETSC_SUCCESS);
  ASSERT_EQ(VecDuplicate(first, &second), PETSC_SUCCESS);
  ASSERT_EQ(VecDuplicate(first, &residual), PETSC_SUCCESS);
  ASSERT_EQ(VecDuplicate(first, &reference), PETSC_SUCCESS);
  ASSERT_EQ(VecSet(first, 1), PETSC_SUCCESS);
  ASSERT_EQ(VecSet(second, 2), PETSC_SUCCESS);
#if PETSC_VERSION_GE(3, 20, 0)
  PetscObjectState firstState = 0, secondState = 0;
  ASSERT_EQ(VecGetState(first, &firstState), PETSC_SUCCESS);
  ASSERT_EQ(VecGetState(second, &secondState), PETSC_SUCCESS);
  ASSERT_EQ(firstState, secondState);
#endif

  auto check = [&](::Vec iterate, size_t expectedUpdates) {
    EXPECT_EQ(SNESComputeFunction(snes.getHandle(), iterate, residual), PETSC_SUCCESS);
    EXPECT_EQ(updates, expectedUpdates);
    // Bypass SNES's caches to obtain the current residual independently.
    state.setData(iterate);
    problem.assemble(AssemblyTarget::RHS);
    EXPECT_EQ(VecCopy(problem.getLinearSystem().getVector(), reference), PETSC_SUCCESS);
    EXPECT_EQ(VecScale(reference, -1), PETSC_SUCCESS);
    PetscBool equal = PETSC_FALSE;
    EXPECT_EQ(VecEqual(residual, reference, &equal), PETSC_SUCCESS);
    EXPECT_EQ(equal, PETSC_TRUE);
#if PETSC_VERSION_GE(3, 20, 0)
    const auto matrix = problem.getLinearSystem().getOperator();
    EXPECT_EQ(
      SNESComputeJacobian(snes.getHandle(), iterate, matrix, matrix), PETSC_SUCCESS);
    EXPECT_EQ(updates, expectedUpdates);
    ::Mat assembled = nullptr;
    EXPECT_EQ(MatDuplicate(matrix, MAT_COPY_VALUES, &assembled), PETSC_SUCCESS);
    problem.assemble(AssemblyTarget::LHS);
    EXPECT_EQ(MatEqual(matrix, assembled, &equal), PETSC_SUCCESS);
    EXPECT_EQ(equal, PETSC_TRUE);
    EXPECT_EQ(MatDestroy(&assembled), PETSC_SUCCESS);
#endif
  };
  check(first, 1);
  check(second, 2);
  check(first, 3);
#if PETSC_VERSION_GE(3, 20, 0)
  check(first, 3); // Same object and state: no new synchronization is needed.
#else
  check(first, 4); // PETSc 3.19 deliberately uses the uncached path.
#endif
  ASSERT_EQ(VecSet(first, 3), PETSC_SUCCESS);
#if PETSC_VERSION_GE(3, 20, 0)
  check(first, 4);
#else
  check(first, 5);
#endif

  EXPECT_EQ(VecDestroy(&reference), PETSC_SUCCESS);
  EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
  EXPECT_EQ(VecDestroy(&second), PETSC_SUCCESS);
  EXPECT_EQ(VecDestroy(&first), PETSC_SUCCESS);
}

/// @brief A manufactured constant root exercises the complete nonlinear solve.
TEST(Rodin_PETSc_SNESState, CubicManufacturedRootAndWarmStart)
{
  auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {3});
  P1<PetscScalar> space(mesh);
  PETSc::Variational::GridFunction state(space);
  state = Zero();
  PETSc::Variational::TrialFunction du(space);
  PETSc::Variational::TestFunction v(space);
  Problem problem(du, v);
#ifdef PETSC_USE_COMPLEX
  const auto one = ComplexFunction(1);
#else
  const auto one = RealFunction(1);
#endif
  problem = Integral((one + 3 * state * state) * du, v) +
    Integral(state + state * state * state - 2 * one, v);
  problem.assemble();
  Solver::KSP ksp(problem);
  Solver::SNES snes(ksp);
  snes.setType(SNESNEWTONLS).setTolerances(1e-12, 1e-12, 1e-12, 30, 200);
  snes.setStateUpdate([&](const PETSc::Math::Vector& x) { state.setData(x); });
  snes.solve();
  EXPECT_TRUE(snes.converged());
  EXPECT_GT(snes.getConvergedReason(), 0);
  EXPECT_GT(snes.getIterationNumber(), 0);
  ::Vec error = nullptr;
  const auto solution = problem.getLinearSystem().getSolution();
  ASSERT_EQ(VecDuplicate(solution, &error), PETSC_SUCCESS);
  ASSERT_EQ(VecCopy(solution, error), PETSC_SUCCESS);
  ASSERT_EQ(VecShift(error, -1), PETSC_SUCCESS);
  PetscReal norm = 0;
  ASSERT_EQ(VecNorm(error, NORM_INFINITY, &norm), PETSC_SUCCESS);
  EXPECT_LT(norm, 1e-10);
  EXPECT_EQ(VecDestroy(&error), PETSC_SUCCESS);
  snes.solve();
  EXPECT_TRUE(snes.converged());
  EXPECT_EQ(snes.getIterationNumber(), 0);
}

#ifdef PETSC_USE_COMPLEX
/**
 * @brief Newton solves recover the non-real manufactured root and preserve warm starts.
 *
 * The residual is @f$F(u)=u+u^3-(z+z^3)@f$, with
 * @f$z=1+i/2@f$. Its tangent is @f$1+3u^2@f$.
 */
TEST(Rodin_PETSc_SNESState, CubicNonrealManufacturedRootAndWarmStart)
{
  constexpr Real residualTolerance = 1e-12;
  constexpr Real solutionTolerance = 1e-10;
  const PetscScalar exact(1, 0.5);
  const PetscScalar initial(0.8, 0.4);
  auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {3});
  P1<PetscScalar> space(mesh);
  PETSc::Variational::GridFunction state(space);
  state = ComplexFunction(initial);
  PETSc::Variational::TrialFunction du(space);
  PETSc::Variational::TestFunction v(space);
  const auto one = ComplexFunction(1);
  const auto forcing = ComplexFunction(exact + exact * exact * exact);
  Problem problem(du, v);
  problem = Integral((one + 3 * state * state) * du, v) +
    Integral(state + state * state * state - forcing, v);
  problem.assemble();
  const auto solution = problem.getLinearSystem().getSolution();
  ASSERT_EQ(VecSet(solution, initial), PETSC_SUCCESS);
  Solver::KSP ksp(problem);
  Solver::SNES snes(ksp);
  snes.setType(SNESNEWTONLS)
    .setTolerances(residualTolerance, residualTolerance, residualTolerance, 30, 200);
  snes.setStateUpdate([&](const PETSc::Math::Vector& x) { state.setData(x); });
  snes.solve();
  EXPECT_TRUE(snes.converged());
  EXPECT_GT(snes.getConvergedReason(), 0);
  EXPECT_GT(snes.getIterationNumber(), 0);
  ::Vec error = nullptr;
  ASSERT_EQ(VecDuplicate(solution, &error), PETSC_SUCCESS);
  ASSERT_EQ(VecCopy(solution, error), PETSC_SUCCESS);
  ASSERT_EQ(VecShift(error, -exact), PETSC_SUCCESS);
  PetscReal norm = 0;
  ASSERT_EQ(VecNorm(error, NORM_INFINITY, &norm), PETSC_SUCCESS);
  EXPECT_LT(norm, solutionTolerance);
  EXPECT_EQ(VecDestroy(&error), PETSC_SUCCESS);
  snes.solve();
  EXPECT_TRUE(snes.converged());
  EXPECT_EQ(snes.getIterationNumber(), 0);
}
#endif

#ifndef PETSC_USE_COMPLEX
/**
 * @brief Domain validation rejects a negative Newton trial before coefficient evaluation.
 *
 * The manufactured equation is @f$1/u-4=0@f$ on @f$u>0@f$.
 * Starting from one gives the inadmissible full Newton step @f$u=-2@f$.
 */
TEST(Rodin_PETSc_SNESState, DomainValidationBacktracksBeforeAssembly)
{
  auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Segment, {2});
  P0<Real> space(mesh);
  PETSc::Variational::GridFunction state(space);
  state = Real(1);
  PETSc::Variational::TrialFunction du(space);
  PETSc::Variational::TestFunction v(space);
  du.getSolution() = state;
  size_t invalidEvaluations = 0, rejected = 0;
  const auto inverse = RealFunction([&](const Point& point) {
    const Real value = state(point);
    if (value <= 0)
      ++invalidEvaluations;
    return 1 / value;
  });
  Problem problem(du, v);
  problem = Integral(-inverse * inverse * du, v) + Integral(inverse - 4, v);
  problem.assemble();
  Solver::KSP ksp(problem);
  Solver::SNES snes(ksp);
  snes.setTolerances(1e-12, 0, 0, 30, 200)
    .setStateUpdate([&](const PETSc::Math::Vector& x) { state.setData(x); })
    .setStateValidation([&] {
      Index index;
      const bool valid = state.min(index) > 0;
      if (!valid)
        ++rejected;
      return valid;
    });
  snes.solve();
  EXPECT_TRUE(snes.converged());
  EXPECT_GT(rejected, 0);
  EXPECT_EQ(invalidEvaluations, 0);
  Index index;
  EXPECT_NEAR(state.min(index), 0.25, 1e-12);
  EXPECT_LT(snes.getResidualNorm(), 1e-12);
}
#endif

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
