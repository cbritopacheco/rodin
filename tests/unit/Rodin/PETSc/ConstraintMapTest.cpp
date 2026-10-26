/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <petsc.h>
#include "Rodin/PETSc/Assembly/ConstraintMap.h"

using namespace Rodin;

/// @brief Remote masters and affine offsets are resolved on ranks with no local rows too.
TEST(PETScConstraintMap, DistributedAffineClosure)
{
  int rank = 0, count = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &count);
  for (const bool fixedTerminal : {false, true})
  {
    PETSc::Assembly::ConstraintMap constraints(3);
    using Entry = PETSc::Assembly::ConstraintMap::Entry;
    if (rank == 0)
    {
      constraints.setIdentification(
        0, std::vector<Entry>{{1, PetscScalar(2)}}, PetscScalar(1));
    }
    if (rank == (count > 1 ? 1 : 0))
    {
      constraints.setIdentification(
        1, std::vector<Entry>{{2, PetscScalar(3)}}, PetscScalar(4));
    }
    if (fixedTerminal && rank == count - 1)
      constraints.setFixed(2, PetscScalar(5));
    constraints.finalize(PETSC_COMM_WORLD);
    if (fixedTerminal)
    {
      ASSERT_TRUE(constraints.isFixed(0));
      ASSERT_TRUE(constraints.isFixed(1));
      EXPECT_EQ(constraints.getFixedValue(0), PetscScalar(39));
      EXPECT_EQ(constraints.getFixedValue(1), PetscScalar(19));
    }
    else
    {
      ASSERT_EQ(constraints.expand(0).size(), 1u);
      EXPECT_EQ(constraints.expand(0)[0].index, 2u);
      EXPECT_EQ(constraints.expand(0)[0].coefficient, PetscScalar(6));
      EXPECT_EQ(constraints.getIdentificationValue(0), PetscScalar(9));
    }
  }
}

/// @brief Every rank diagnoses a cycle whose edges are supplied by different ranks.
TEST(PETScConstraintMap, DistributedCycle)
{
  int rank = 0, count = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &count);
  PETSc::Assembly::ConstraintMap constraints(2);
  using Entry = PETSc::Assembly::ConstraintMap::Entry;
  if (rank == 0)
    constraints.setIdentification(0, std::vector<Entry>{{1, PetscScalar(1)}});
  if (rank == (count > 1 ? 1 : 0))
    constraints.setIdentification(1, std::vector<Entry>{{0, PetscScalar(1)}});
  EXPECT_THROW(constraints.finalize(PETSC_COMM_WORLD), Alert::Exception);
}

/// @brief Equivalent replicated rows do not become duplicate local identifications.
TEST(PETScConstraintMap, ReplicatedRows)
{
  PETSc::Assembly::ConstraintMap constraints(3);
  using Entry = PETSc::Assembly::ConstraintMap::Entry;
  constraints.setIdentification(
    0, std::vector<Entry>{{1, PetscScalar(2)}, {2, PetscScalar(3)}});
  constraints.setFixed(1, PetscScalar(4));
  constraints.finalize(PETSC_COMM_WORLD);
  ASSERT_EQ(constraints.expand(0).size(), 1u);
  EXPECT_EQ(constraints.expand(0)[0].index, 2u);
  EXPECT_EQ(constraints.expand(0)[0].coefficient, PetscScalar(3));
  EXPECT_EQ(constraints.getIdentificationValue(0), PetscScalar(8));
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
