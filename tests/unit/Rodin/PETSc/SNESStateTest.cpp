/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief SNES iterate identity and state-synchronization contract. */

#include <gtest/gtest.h>
#include <petsc.h>

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

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
