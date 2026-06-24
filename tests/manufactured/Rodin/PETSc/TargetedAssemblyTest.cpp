/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>

#include <gtest/gtest.h>

#include <petsc.h>

#include <Rodin/Assembly.h>
#include <Rodin/Configure.h>
#include <Rodin/Geometry.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <template <class, class> class Assembler>
  PETSc::Math::LinearSystem assembleLocal(
      Variational::AssemblyTarget target,
      bool targeted)
  {
    auto mesh =
      Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType = Problem<LinearSystemType, decltype(u), decltype(v)>;
    typename ProblemType::ProblemBodyType body =
      Integral(u, v) - Integral(RealFunction(1.0), v);

    Assembly::ProblemAssemblyInput<
      typename ProblemType::ProblemBodyType,
      decltype(u), decltype(v)> input(body, u, v);

    LinearSystemType ls(PETSC_COMM_SELF);
    Assembler<LinearSystemType, ProblemType> assembler;
    if (targeted)
      assembler.execute(ls, input, target);
    else
      assembler.execute(ls, input);
    return ls;
  }

  void expectSameVector(Vec expected, Vec actual)
  {
    PetscErrorCode ierr;
    PetscInt nExpected = 0;
    PetscInt nActual = 0;
    ierr = VecGetSize(expected, &nExpected);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetSize(actual, &nActual);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ASSERT_EQ(nExpected, nActual);

    for (PetscInt i = 0; i < nExpected; i++)
    {
      PetscScalar a = 0;
      PetscScalar b = 0;
      ierr = VecGetValues(expected, 1, &i, &a);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecGetValues(actual, 1, &i, &b);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_LE(PetscAbsScalar(a - b), 1e-14) << "row " << i;
    }
  }

  void expectSameMatrix(Mat expected, Mat actual)
  {
    PetscErrorCode ierr;
    PetscInt expectedRows = 0;
    PetscInt expectedCols = 0;
    PetscInt actualRows = 0;
    PetscInt actualCols = 0;
    ierr = MatGetSize(expected, &expectedRows, &expectedCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatGetSize(actual, &actualRows, &actualCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ASSERT_EQ(expectedRows, actualRows);
    ASSERT_EQ(expectedCols, actualCols);

    for (PetscInt i = 0; i < expectedRows; i++)
    {
      for (PetscInt j = 0; j < expectedCols; j++)
      {
        PetscScalar a = 0;
        PetscScalar b = 0;
        ierr = MatGetValues(expected, 1, &i, 1, &j, &a);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        ierr = MatGetValues(actual, 1, &i, 1, &j, &b);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        EXPECT_LE(PetscAbsScalar(a - b), 1e-14)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

  template <template <class, class> class Assembler>
  void checkLocalBackendTargetedAssembly()
  {
    auto full = assembleLocal<Assembler>(
        Variational::AssemblyTarget::LHS, false);
    auto lhs = assembleLocal<Assembler>(
        Variational::AssemblyTarget::LHS, true);
    auto rhs = assembleLocal<Assembler>(
        Variational::AssemblyTarget::RHS, true);

    expectSameMatrix(full.getOperator(), lhs.getOperator());
    expectSameVector(full.getVector(), rhs.getVector());
  }

  TEST(PETSc_TargetedAssembly, SequentialBackendAssemblesLHSAndRHS)
  {
    checkLocalBackendTargetedAssembly<Assembly::Sequential>();
  }

#ifdef RODIN_USE_OPENMP
  TEST(PETSc_TargetedAssembly, OpenMPBackendAssemblesLHSAndRHS)
  {
    checkLocalBackendTargetedAssembly<Assembly::OpenMP>();
  }
#endif
}

int main(int argc, char** argv)
{
  [[maybe_unused]] PetscErrorCode ierr =
    PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
