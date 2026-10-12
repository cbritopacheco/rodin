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
  /// @brief Verifies eliminate does not write locked solution vector for PET sc linear system by checking tolerance-based numerical results, exact expected values.
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

/// @brief Raw PETSc systems share referenced handles on copy and transfer them on move.
TEST(PETSc_LinearSystem, RawHandleCopyMoveAndAssignmentLifetime)
{
  using System = Rodin::Math::LinearSystem<::Mat, ::Vec>;
  ::Mat matrix = nullptr;
  ::Vec solution = nullptr;
  System survivor(PETSC_COMM_SELF);
  {
    System original(PETSC_COMM_SELF);
    matrix = original.getOperator();
    solution = original.getSolution();
    ASSERT_EQ(MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, 2, 2), PETSC_SUCCESS);
    ASSERT_EQ(MatSetType(matrix, MATSEQAIJ), PETSC_SUCCESS);
    ASSERT_EQ(MatSetUp(matrix), PETSC_SUCCESS);
    ASSERT_EQ(VecSetSizes(solution, PETSC_DECIDE, 2), PETSC_SUCCESS);
    ASSERT_EQ(VecSetType(solution, VECSEQ), PETSC_SUCCESS);
    ASSERT_EQ(VecSet(solution, 7), PETSC_SUCCESS);
    System copied(original);
    EXPECT_EQ(copied.getOperator(), matrix);
    EXPECT_EQ(copied.getSolution(), solution);
    EXPECT_EQ(copied.getVector(), original.getVector());
    System moved(std::move(copied));
    EXPECT_EQ(moved.getOperator(), matrix);
    EXPECT_EQ(copied.getOperator(), nullptr);
    EXPECT_EQ(copied.getSolution(), nullptr);
    EXPECT_EQ(copied.getVector(), nullptr);
    survivor = moved;
    auto& same = survivor;
    survivor = same;
    EXPECT_EQ(survivor.getOperator(), matrix);
    System assigned(PETSC_COMM_SELF);
    assigned = std::move(moved);
    EXPECT_EQ(moved.getOperator(), nullptr);
    EXPECT_EQ(assigned.getOperator(), matrix);
    survivor = std::move(assigned);
    EXPECT_EQ(assigned.getSolution(), nullptr);
  }
  EXPECT_EQ(survivor.getOperator(), matrix);
  EXPECT_EQ(survivor.getSolution(), solution);
  PetscReal norm = 0;
  ASSERT_EQ(VecNorm(survivor.getSolution(), NORM_INFINITY, &norm), PETSC_SUCCESS);
  EXPECT_DOUBLE_EQ(norm, 7);
  System movedAgain(std::move(survivor));
  System copiedFromEmpty(survivor);
  EXPECT_EQ(copiedFromEmpty.getOperator(), nullptr);
  EXPECT_EQ(copiedFromEmpty.getSolution(), nullptr);
  EXPECT_EQ(copiedFromEmpty.getVector(), nullptr);
}
