/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <array>
#include <string>
#include <vector>

#include <petsc.h>

#include <Rodin/Types.h>
#include <Rodin/Alert/Exception.h>
#include <Rodin/PETSc/Math/LinearSystem.h>

namespace
{
  /// @brief Unsupported periodic merging raises a contextual Alert without replacing handles.
  TEST(PETSc_LinearSystem, UnsupportedMergeReportsMemberContext)
  {
    Rodin::PETSc::Math::LinearSystem system(PETSC_COMM_SELF);
    const auto matrix = system.getOperator();
    const auto vector = system.getVector();
    const auto solution = system.getSolution();
    Rodin::IndexMap<std::pair<Rodin::IndexArray, Rodin::Math::Vector<Rodin::Real>>> dofs;
    try
    {
      system.merge(dofs);
      FAIL() << "Unsupported periodic merging must be rejected.";
    }
    catch (const Rodin::Alert::Exception& error)
    {
      EXPECT_NE(std::string(error.what()).find("merge"), std::string::npos);
    }
    EXPECT_EQ(system.getOperator(), matrix);
    EXPECT_EQ(system.getVector(), vector);
    EXPECT_EQ(system.getSolution(), solution);
  }

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

/// @brief Field-split copies retain independent names and index sets after source changes.
TEST(PETSc_LinearSystem, FieldSplitsCopyRetainsIndependentIndices)
{
  using Fields = Rodin::PETSc::Math::LinearSystem::FieldSplits;
  const std::array<PetscInt, 3> velocity{0, 2, 4};
  const std::array<PetscInt, 2> pressure{1, 3};
  ::IS velocitySet = nullptr, pressureSet = nullptr;
  ASSERT_EQ(ISCreateGeneral(PETSC_COMM_SELF, velocity.size(), velocity.data(),
              PETSC_COPY_VALUES, &velocitySet),
    PETSC_SUCCESS);
  ASSERT_EQ(ISCreateGeneral(PETSC_COMM_SELF, pressure.size(), pressure.data(),
              PETSC_COPY_VALUES, &pressureSet),
    PETSC_SUCCESS);
  Fields original({{"velocity", velocitySet}, {"pressure", pressureSet}});
  const Fields copied(original);
  ASSERT_EQ(copied.getSplits().size(), 2u);
  EXPECT_NE(copied.getSplits()[0].is, velocitySet);
  EXPECT_NE(copied.getSplits()[1].is, pressureSet);
  const std::array<PetscInt, 2> replacement{2, 4};
  ASSERT_EQ(ISGeneralSetIndices(original.getSplits()[0].is, replacement.size(),
              replacement.data(), PETSC_COPY_VALUES),
    PETSC_SUCCESS);
  original.getSplits()[0].name = "changed";
  original.clear();
  EXPECT_TRUE(original.getSplits().empty());
  const std::array<std::vector<PetscInt>, 2> expected{
    std::vector<PetscInt>{0, 2, 4}, std::vector<PetscInt>{1, 3}};
  const std::array<std::string, 2> names{"velocity", "pressure"};
  for (size_t i = 0; i < expected.size(); ++i)
  {
    const auto& split = copied.getSplits()[i];
    EXPECT_EQ(split.name, names[i]);
    PetscInt size = 0;
    ASSERT_EQ(ISGetLocalSize(split.is, &size), PETSC_SUCCESS);
    ASSERT_EQ(size, expected[i].size());
    const PetscInt* indices = nullptr;
    ASSERT_EQ(ISGetIndices(split.is, &indices), PETSC_SUCCESS);
    for (size_t j = 0; j < expected[i].size(); ++j)
      EXPECT_EQ(indices[j], expected[i][j]);
    EXPECT_EQ(ISRestoreIndices(split.is, &indices), PETSC_SUCCESS);
  }
}

/// @brief Copy assignment releases previous owned handles and is safe under self-assignment.
TEST(PETSc_LinearSystem, FieldSplitAssignmentReleasesAndReplacesOwnedSets)
{
  using Fields = Rodin::PETSc::Math::LinearSystem::FieldSplits;
  const std::array<PetscInt, 2> indices{1, 4};
  ::IS source = nullptr, observer = nullptr;
  ASSERT_EQ(ISCreateGeneral(PETSC_COMM_SELF, indices.size(), indices.data(),
              PETSC_COPY_VALUES, &source),
    PETSC_SUCCESS);
  ASSERT_EQ(ISCreateStride(PETSC_COMM_SELF, 3, 0, 2, &observer), PETSC_SUCCESS);
  Fields original({{"selected", source}});
  Fields assigned({{"previous", observer}});
  ASSERT_EQ(PetscObjectReference(reinterpret_cast<PetscObject>(observer)), PETSC_SUCCESS);
  EXPECT_EQ(&(assigned = original), &assigned);
  PetscInt references = 0;
  ASSERT_EQ(PetscObjectGetReference(reinterpret_cast<PetscObject>(observer), &references),
    PETSC_SUCCESS);
  EXPECT_EQ(references, 1);
  ASSERT_EQ(ISDestroy(&observer), PETSC_SUCCESS);
  ASSERT_EQ(assigned.getSplits().size(), 1u);
  EXPECT_NE(assigned.getSplits()[0].is, source);
  const auto handle = assigned.getSplits()[0].is;
  auto& same = assigned;
  EXPECT_EQ(&(assigned = same), &assigned);
  EXPECT_EQ(assigned.getSplits()[0].is, handle);
  original.clear();
  const Fields& view = assigned;
  EXPECT_EQ(view.getSplits()[0].name, "selected");
  PetscInt count = 0;
  ASSERT_EQ(ISGetLocalSize(handle, &count), PETSC_SUCCESS);
  ASSERT_EQ(count, indices.size());
  const PetscInt* values = nullptr;
  ASSERT_EQ(ISGetIndices(handle, &values), PETSC_SUCCESS);
  for (size_t i = 0; i < indices.size(); ++i)
    EXPECT_EQ(values[i], indices[i]);
  EXPECT_EQ(ISRestoreIndices(handle, &values), PETSC_SUCCESS);
  assigned.clear();
  EXPECT_TRUE(assigned.getSplits().empty());
}

/// @brief Public default descriptors retain their unbound value through copying.
TEST(PETSc_LinearSystem, DefaultFieldSplitDescriptorsCopyWithoutHandles)
{
  using Fields = Rodin::PETSc::Math::LinearSystem::FieldSplits;
  Fields empty;
  const Fields emptyCopy(empty);
  EXPECT_TRUE(emptyCopy.getSplits().empty());
  Fields original({Fields::Split{}});
  const Fields copied(original);
  Fields assigned;
  assigned = original;
  original.clear();
  for (const Fields* fields : {&copied, static_cast<const Fields*>(&assigned)})
  {
    ASSERT_EQ(fields->getSplits().size(), 1u);
    EXPECT_TRUE(fields->getSplits()[0].name.empty());
    EXPECT_EQ(fields->getSplits()[0].is, PETSC_NULLPTR);
  }
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
