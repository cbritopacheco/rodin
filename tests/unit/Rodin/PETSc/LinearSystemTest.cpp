/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <array>

#include <petsc.h>

#include <Rodin/Types.h>
#include <Rodin/PETSc/Math/LinearSystem.h>

namespace
{
  TEST(PETSc_LinearSystem, EliminateDoesNotWriteLockedSolutionVector)
  {
    Rodin::PETSc::Math::LinearSystem system(PETSC_COMM_SELF);

    auto& A = system.getOperator();
    auto& b = system.getVector();
    auto& x = system.getSolution();

    PetscErrorCode ierr;
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 3, 3);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetFromOptions(A);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetUp(A);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    ierr = VecSetSizes(b, PETSC_DECIDE, 3);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(b);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    ierr = VecSetSizes(x, PETSC_DECIDE, 3);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(x);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    for (PetscInt i = 0; i < 3; ++i)
    {
      for (PetscInt j = 0; j < 3; ++j)
      {
        const PetscScalar value = (i == j) ? 2.0 : 0.25;
        ierr = MatSetValue(A, i, j, value, INSERT_VALUES);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
      }
      ierr = VecSetValue(b, i, 1.0, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecSetValue(x, i, 10.0 + static_cast<PetscScalar>(i), INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyBegin(b);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(b);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyBegin(x);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(x);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    Rodin::IndexMap<PetscScalar> dofs;
    dofs.emplace(1, 42.0);

    ierr = VecLockReadPush(x);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    system.eliminate(dofs);
    ierr = VecLockReadPop(x);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    PetscScalar value = 0;
    ierr = VecGetValues(x, 1, std::array<PetscInt, 1>{1}.data(), &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(value)), 11.0);

    PetscInt row = 1;
    PetscInt col = 1;
    ierr = MatGetValues(A, 1, &row, 1, &col, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(value)), 1.0);

    ierr = VecGetValues(b, 1, &row, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(value)), 42.0);
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
