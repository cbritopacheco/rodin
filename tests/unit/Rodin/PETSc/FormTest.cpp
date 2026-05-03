/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <mpi.h>
#include <petsc.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  TEST(PETSc_Form, SequentialLinearFormUsesSelfCommunicatorAndAssembles)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();

    auto& b = lf.getVector();

    MPI_Comm comm = MPI_COMM_NULL;
    PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(b), &comm);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    int commSize = 0;
    MPI_Comm_size(comm, &commSize);
    EXPECT_EQ(commSize, 1);

    PetscInt localSize = 0;
    PetscInt globalSize = 0;
    ierr = VecGetLocalSize(b, &localSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetSize(b, &globalSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(localSize, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalSize, static_cast<PetscInt>(fes.getSize()));
  }

  TEST(PETSc_Form, SequentialBilinearFormUsesSelfCommunicatorAndAssembles)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    auto& A = bf.getOperator();

    MPI_Comm comm = MPI_COMM_NULL;
    PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(A), &comm);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    int commSize = 0;
    MPI_Comm_size(comm, &commSize);
    EXPECT_EQ(commSize, 1);

    PetscInt localRows = 0;
    PetscInt localCols = 0;
    PetscInt globalRows = 0;
    PetscInt globalCols = 0;
    ierr = MatGetLocalSize(A, &localRows, &localCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatGetSize(A, &globalRows, &globalCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(localRows, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(localCols, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalRows, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalCols, static_cast<PetscInt>(fes.getSize()));
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
