/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Targeted assembly test manufactured regression tests.
 *
 * These tests assemble Rodin variational forms for a Targeted Assembly Test manufactured regression, solve the problem on the configured mesh, and compare against analytic fields or expected residual/error behavior. They protect the PETSc-backed assembly and solve path, including boundary-condition handling, geometry coverage, and numerical accuracy of the manufactured workflow.
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
    Variational::AssemblyTarget target, bool targeted)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType = Problem<LinearSystemType, decltype(u), decltype(v)>;
    typename ProblemType::ProblemBodyType body =
      Integral(u, v) - Integral(RealFunction(1.0), v);

    Assembly::ProblemAssemblyInput<typename ProblemType::ProblemBodyType, decltype(u),
      decltype(v)>
      input(body, u, v);

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

  void expectVectorConstant(Vec vector, PetscScalar expected)
  {
    PetscErrorCode ierr;
    PetscInt n = 0;
    ierr = VecGetSize(vector, &n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    for (PetscInt i = 0; i < n; i++)
    {
      PetscScalar value = 0;
      ierr = VecGetValues(vector, 1, &i, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_LE(PetscAbsScalar(value - expected), 1e-14) << "row " << i;
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
        EXPECT_LE(PetscAbsScalar(a - b), 1e-14) << "entry (" << i << ", " << j << ")";
      }
    }
  }

  template <template <class, class> class Assembler>
  void checkLocalBackendTargetedAssembly()
  {
    auto full = assembleLocal<Assembler>(Variational::AssemblyTarget::LHS, false);
    auto lhs = assembleLocal<Assembler>(Variational::AssemblyTarget::LHS, true);
    auto rhs = assembleLocal<Assembler>(Variational::AssemblyTarget::RHS, true);

    expectSameMatrix(full.getOperator(), lhs.getOperator());
    expectSameVector(full.getVector(), rhs.getVector());
  }

  // Number of stored nonzeros, i.e. the size of the PETSc-owned nonzero
  // pattern under the current PETSc options.
  PetscInt matrixNonzeros(Mat A)
  {
    MatInfo info;
    PetscErrorCode ierr = MatGetInfo(A, MAT_LOCAL, &info);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    return static_cast<PetscInt>(info.nz_used);
  }

  // Cumulative number of mallocs performed during MatSetValues over the matrix
  // lifetime. A reused pattern triggers no additional mallocs.
  PetscReal matrixMallocs(Mat A)
  {
    MatInfo info;
    PetscErrorCode ierr = MatGetInfo(A, MAT_LOCAL, &info);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    return info.mallocs;
  }

  // Assemble a P1 stiffness system (Grad . Grad scaled by gammaValue) on an
  // n-by-n triangle grid into the supplied linear system. Re-running with the
  // same n exercises the MatrixSetup reuse path; changing n forces a full
  // structural setup.
  template <template <class, class> class Assembler>
  void assembleStiffness(PETSc::Math::LinearSystem& ls, size_t n, double gammaValue,
    bool emptyLHS, double sourceValue)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    RealFunction gamma(gammaValue);

    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType = Problem<LinearSystemType, decltype(u), decltype(v)>;
    auto lhs = Integral(gamma * Grad(u), Grad(v));
    if (emptyLHS)
      lhs.over(999);

    typename ProblemType::ProblemBodyType body =
      lhs - Integral(RealFunction(sourceValue), v);

    Assembly::ProblemAssemblyInput<typename ProblemType::ProblemBodyType, decltype(u),
      decltype(v)>
      input(body, u, v);

    Assembler<LinearSystemType, ProblemType> assembler;
    assembler.execute(ls, input);
  }

  template <template <class, class> class Assembler>
  void assembleStiffness(
    PETSc::Math::LinearSystem& ls, size_t n, double gammaValue, bool emptyLHS)
  {
    assembleStiffness<Assembler>(ls, n, gammaValue, emptyLHS, 1.0);
  }

  template <template <class, class> class Assembler>
  void assembleStiffness(PETSc::Math::LinearSystem& ls, size_t n, double gammaValue)
  {
    assembleStiffness<Assembler>(ls, n, gammaValue, false, 1.0);
  }

  // Re-assembling the identical problem into the same matrix must reuse the
  // nonzero pattern: identical dimensions, identical nonzero count and no extra
  // mallocs, so PETSc reports SAME_NONZERO_PATTERN.
  template <template <class, class> class Assembler>
  void checkReassemblyKeepsPattern()
  {
    PETSc::Math::LinearSystem ls(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(ls, 4, 1.0);

    PetscInt rows1 = 0;
    PetscInt cols1 = 0;
    PetscErrorCode ierr = MatGetSize(ls.getOperator(), &rows1, &cols1);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    const PetscInt nz1 = matrixNonzeros(ls.getOperator());
    const PetscReal mallocs1 = matrixMallocs(ls.getOperator());

    assembleStiffness<Assembler>(ls, 4, 1.0);

    PetscInt rows2 = 0;
    PetscInt cols2 = 0;
    ierr = MatGetSize(ls.getOperator(), &rows2, &cols2);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(rows1, rows2);
    EXPECT_EQ(cols1, cols2);
    EXPECT_EQ(nz1, matrixNonzeros(ls.getOperator()));
    EXPECT_EQ(mallocs1, matrixMallocs(ls.getOperator()));
  }

  // A problem with different global dimensions must produce a different nonzero
  // pattern. Each distinct-size problem uses its own matrix object (the
  // production path for adaptive remeshing); the two patterns must differ in
  // both dimension and stored nonzero count.
  //
  // Note: re-assembling a *resized* problem into the *same* matrix object is
  // intentionally not exercised here. A LinearSystem is tied to fixed finite
  // element spaces; a different space or mesh must use a fresh LinearSystem.
  template <template <class, class> class Assembler>
  void checkReassemblyChangesPattern()
  {
    PETSc::Math::LinearSystem coarse(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(coarse, 4, 1.0);

    PetscInt rows1 = 0;
    PetscInt cols1 = 0;
    PetscErrorCode ierr = MatGetSize(coarse.getOperator(), &rows1, &cols1);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    const PetscInt nz1 = matrixNonzeros(coarse.getOperator());

    PETSc::Math::LinearSystem fine(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(fine, 6, 1.0);

    PetscInt rows2 = 0;
    PetscInt cols2 = 0;
    ierr = MatGetSize(fine.getOperator(), &rows2, &cols2);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_NE(rows1, rows2);
    EXPECT_NE(cols1, cols2);
    EXPECT_NE(nz1, matrixNonzeros(fine.getOperator()));
  }

  // The right-hand side vector is rebuilt from scratch on every assembly and
  // must be zeroed on reuse. The solution vector is a solver iterate and must
  // not be zeroed when the same LinearSystem is assembled again.
  template <template <class, class> class Assembler>
  void checkVectorReuse()
  {
    PETSc::Math::LinearSystem ls(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(ls, 4, 1.0, false, 1.0);

    PetscErrorCode ierr = VecSet(ls.getSolution(), 3.0);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    assembleStiffness<Assembler>(ls, 4, 1.0, false, 2.0);
    expectVectorConstant(ls.getSolution(), 3.0);

    PETSc::Math::LinearSystem expected(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(expected, 4, 1.0, false, 2.0);
    expectSameVector(expected.getVector(), ls.getVector());
  }

  // The PETSc policy may also be changed between two assemblies of the same
  // LinearSystem. If the first assembly leaves an empty matrix pattern, then
  // explicitly allowing new nonzero allocations before reassembly must let
  // PETSc grow the pattern instead of being overridden by MatrixSetup.
  template <template <class, class> class Assembler>
  void checkNewNonzeroOptionHonoredBeforeReassembly()
  {
    PETSc::Math::LinearSystem ls(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(ls, 4, 1.0, true);
    ASSERT_EQ(matrixNonzeros(ls.getOperator()), 0);

    PetscErrorCode ierr =
      MatSetOption(ls.getOperator(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    assembleStiffness<Assembler>(ls, 4, 1.0);
    EXPECT_GT(matrixNonzeros(ls.getOperator()), 0);
  }

#ifndef NDEBUG
  // Symmetrically, explicitly forbidding new nonzero allocations before the
  // same reassembly must be enforced by PETSc. In debug builds the Rodin PETSc
  // backend checks PETSc errors with assert, so a hidden override to allow new
  // nonzeros would make this death test fail.
  template <template <class, class> class Assembler>
  void checkNewNonzeroOptionForbidsPatternGrowthBeforeReassembly()
  {
    PETSc::Math::LinearSystem ls(PETSC_COMM_SELF);
    assembleStiffness<Assembler>(ls, 4, 1.0, true);
    ASSERT_EQ(matrixNonzeros(ls.getOperator()), 0);

    PetscErrorCode ierr =
      MatSetOption(ls.getOperator(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_DEATH(assembleStiffness<Assembler>(ls, 4, 1.0), "");
  }
#endif

  /// @brief Verifies sequential backend assembles LHS and RHS for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, SequentialBackendAssemblesLHSAndRHS)
  {
    checkLocalBackendTargetedAssembly<Assembly::Sequential>();
  }

  /// @brief Verifies sequential reassembly keeps nonzero pattern for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, SequentialReassemblyKeepsNonzeroPattern)
  {
    checkReassemblyKeepsPattern<Assembly::Sequential>();
  }

  /// @brief Verifies sequential reassembly changes nonzero pattern for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, SequentialReassemblyChangesNonzeroPattern)
  {
    checkReassemblyChangesPattern<Assembly::Sequential>();
  }

  /// @brief Verifies sequential honors new nonzero option before reassembly for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, SequentialHonorsNewNonzeroOptionBeforeReassembly)
  {
    checkNewNonzeroOptionHonoredBeforeReassembly<Assembly::Sequential>();
  }

  /// @brief Verifies sequential reuses vectors correctly for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, SequentialReusesVectorsCorrectly)
  {
    checkVectorReuse<Assembly::Sequential>();
  }

#ifndef NDEBUG
  /// @brief Verifies sequential forbids new nonzero before reassembly for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, SequentialForbidsNewNonzeroBeforeReassembly)
  {
    checkNewNonzeroOptionForbidsPatternGrowthBeforeReassembly<Assembly::Sequential>();
  }
#endif

#ifdef RODIN_USE_OPENMP
  /// @brief Verifies open MP backend assembles LHS and RHS for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, OpenMPBackendAssemblesLHSAndRHS)
  {
    checkLocalBackendTargetedAssembly<Assembly::OpenMP>();
  }

  /// @brief Verifies open MP reassembly keeps nonzero pattern for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, OpenMPReassemblyKeepsNonzeroPattern)
  {
    checkReassemblyKeepsPattern<Assembly::OpenMP>();
  }

  /// @brief Verifies open MP reassembly changes nonzero pattern for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, OpenMPReassemblyChangesNonzeroPattern)
  {
    checkReassemblyChangesPattern<Assembly::OpenMP>();
  }

  /// @brief Verifies open MP honors new nonzero option before reassembly for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, OpenMPHonorsNewNonzeroOptionBeforeReassembly)
  {
    checkNewNonzeroOptionHonoredBeforeReassembly<Assembly::OpenMP>();
  }

  /// @brief Verifies open MP reuses vectors correctly for PET sc targeted assembly by checking form assembly.
  TEST(PETSc_TargetedAssembly, OpenMPReusesVectorsCorrectly)
  {
    checkVectorReuse<Assembly::OpenMP>();
  }

#endif
}

/// @brief Initializes the test runtime and runs the GoogleTest suite.
int main(int argc, char** argv)
{
  [[maybe_unused]] PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
