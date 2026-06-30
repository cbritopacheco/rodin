/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>

#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static boost::mpi::environment*  g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  Mesh<Context::Local> makeShardableMesh(size_t n = 5)
  {
    auto mesh =
      Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { n, n });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx, size_t n = 5)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);

    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(n);
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }

    return sharder.gather(0);
  }

  // Per-rank stored nonzeros, i.e. the size of this rank's PETSc-owned
  // nonzero pattern under the current PETSc options.
  PetscInt matrixLocalNonzeros(Mat A)
  {
    MatInfo info;
    PetscErrorCode ierr = MatGetInfo(A, MAT_LOCAL, &info);
    EXPECT_EQ(ierr, PETSC_SUCCESS);
    return static_cast<PetscInt>(info.nz_used);
  }

  // Per-rank cumulative mallocs during MatSetValues. A reused pattern triggers
  // no additional mallocs.
  PetscReal matrixLocalMallocs(Mat A)
  {
    MatInfo info;
    PetscErrorCode ierr = MatGetInfo(A, MAT_LOCAL, &info);
    EXPECT_EQ(ierr, PETSC_SUCCESS);
    return info.mallocs;
  }

  // Total stored nonzeros across all ranks (collective).
  PetscInt matrixGlobalNonzeros(Mat A)
  {
    MatInfo info;
    PetscErrorCode ierr = MatGetInfo(A, MAT_GLOBAL_SUM, &info);
    EXPECT_EQ(ierr, PETSC_SUCCESS);
    return static_cast<PetscInt>(info.nz_used);
  }

  void expectSameOwnedVector(Vec expected, Vec actual)
  {
    PetscErrorCode ierr;
    PetscInt begin = 0;
    PetscInt end = 0;
    ierr = VecGetOwnershipRange(expected, &begin, &end);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    PetscInt actualBegin = 0;
    PetscInt actualEnd = 0;
    ierr = VecGetOwnershipRange(actual, &actualBegin, &actualEnd);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ASSERT_EQ(begin, actualBegin);
    ASSERT_EQ(end, actualEnd);

    for (PetscInt i = begin; i < end; i++)
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

  void expectOwnedVectorConstant(Vec vector, PetscScalar expected)
  {
    PetscErrorCode ierr;
    PetscInt begin = 0;
    PetscInt end = 0;
    ierr = VecGetOwnershipRange(vector, &begin, &end);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    for (PetscInt i = begin; i < end; i++)
    {
      PetscScalar value = 0;
      ierr = VecGetValues(vector, 1, &i, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_LE(PetscAbsScalar(value - expected), 1e-14) << "row " << i;
    }
  }

  void expectSameOwnedMatrix(Mat expected, Mat actual)
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

    PetscInt begin = 0;
    PetscInt end = 0;
    ierr = MatGetOwnershipRange(expected, &begin, &end);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    PetscInt actualBegin = 0;
    PetscInt actualEnd = 0;
    ierr = MatGetOwnershipRange(actual, &actualBegin, &actualEnd);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ASSERT_EQ(begin, actualBegin);
    ASSERT_EQ(end, actualEnd);

    for (PetscInt i = begin; i < end; i++)
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
}

namespace Rodin::Tests::Manufactured::PETSc::MPI
{
  namespace PETSc = ::Rodin::PETSc;

  TEST(PETSc_MPI_TargetedAssembly, AssemblesLHSAndRHS)
  {
    const auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fullFES(mesh);
    PETSc::Variational::TrialFunction uFull(fullFES);
    PETSc::Variational::TestFunction  vFull(fullFES);
    Problem full(uFull, vFull);
    full = Integral(uFull, vFull) - Integral(RealFunction(1.0), vFull);
    full.assemble();

    P1<Real, Mesh<Context::MPI>> lhsFES(mesh);
    PETSc::Variational::TrialFunction uLHS(lhsFES);
    PETSc::Variational::TestFunction  vLHS(lhsFES);
    Problem lhs(uLHS, vLHS);
    lhs = Integral(uLHS, vLHS) - Integral(RealFunction(1.0), vLHS);
    lhs.assemble(Variational::AssemblyTarget::LHS);

    P1<Real, Mesh<Context::MPI>> rhsFES(mesh);
    PETSc::Variational::TrialFunction uRHS(rhsFES);
    PETSc::Variational::TestFunction  vRHS(rhsFES);
    Problem rhs(uRHS, vRHS);
    rhs = Integral(uRHS, vRHS) - Integral(RealFunction(1.0), vRHS);
    rhs.assemble(Variational::AssemblyTarget::RHS);

    expectSameOwnedMatrix(
        full.getLinearSystem().getOperator(),
        lhs.getLinearSystem().getOperator());
    expectSameOwnedVector(
        full.getLinearSystem().getVector(),
        rhs.getLinearSystem().getVector());
  }

  // Re-assembling the identical problem into the same distributed matrix must
  // reuse the nonzero pattern on every rank: identical dimensions, identical
  // per-rank nonzero count and no extra mallocs (SAME_NONZERO_PATTERN).
  TEST(PETSc_MPI_TargetedAssembly, ReassemblyKeepsNonzeroPattern)
  {
    const auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);
    RealFunction gamma(1.0);
    Problem p(u, v);
    p = Integral(gamma * Grad(u), Grad(v)) - Integral(RealFunction(1.0), v);

    p.assemble();
    Mat A = p.getLinearSystem().getOperator();
    PetscInt rows1 = 0;
    PetscInt cols1 = 0;
    ASSERT_EQ(MatGetSize(A, &rows1, &cols1), PETSC_SUCCESS);
    const PetscInt nz1 = matrixLocalNonzeros(A);
    const PetscReal mallocs1 = matrixLocalMallocs(A);

    p.assemble();
    A = p.getLinearSystem().getOperator();
    PetscInt rows2 = 0;
    PetscInt cols2 = 0;
    ASSERT_EQ(MatGetSize(A, &rows2, &cols2), PETSC_SUCCESS);

    EXPECT_EQ(rows1, rows2);
    EXPECT_EQ(cols1, cols2);
    EXPECT_EQ(nz1, matrixLocalNonzeros(A));
    EXPECT_EQ(mallocs1, matrixLocalMallocs(A));
  }

  // The PETSc policy may also be changed between two assemblies of the same
  // distributed LinearSystem. If the first assembly leaves an empty matrix
  // pattern, then explicitly allowing new nonzero allocations before
  // reassembly must let PETSc grow the pattern instead of being overridden by
  // the MPI MatrixSetup path.
  TEST(PETSc_MPI_TargetedAssembly, HonorsNewNonzeroOptionBeforeReassembly)
  {
    const auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);
    RealFunction gamma(1.0);

    auto emptyLHS = Integral(gamma * Grad(u), Grad(v));
    emptyLHS.over(999);

    Problem p(u, v);
    p = emptyLHS - Integral(RealFunction(1.0), v);
    p.assemble();

    Mat A = p.getLinearSystem().getOperator();
    ASSERT_EQ(matrixGlobalNonzeros(A), 0);

    ASSERT_EQ(
        MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
        PETSC_SUCCESS);
    p = Integral(gamma * Grad(u), Grad(v)) - Integral(RealFunction(1.0), v);
    p.assemble();

    EXPECT_GT(matrixGlobalNonzeros(A), 0);
  }

  // The distributed right-hand side vector is rebuilt from scratch on every
  // assembly and must be zeroed on reuse. The distributed solution vector is a
  // solver iterate and must not be zeroed when the same Problem is assembled
  // again.
  TEST(PETSc_MPI_TargetedAssembly, ReusesVectorsCorrectly)
  {
    const auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);
    RealFunction gamma(1.0);

    Problem p(u, v);
    p = Integral(gamma * Grad(u), Grad(v))
      - Integral(RealFunction(1.0), v);
    p.assemble();

    Vec x = p.getLinearSystem().getSolution();
    ASSERT_EQ(VecSet(x, 3.0), PETSC_SUCCESS);

    p = Integral(gamma * Grad(u), Grad(v))
      - Integral(RealFunction(2.0), v);
    p.assemble();
    expectOwnedVectorConstant(p.getLinearSystem().getSolution(), 3.0);

    P1<Real, Mesh<Context::MPI>> expectedFES(mesh);
    PETSc::Variational::TrialFunction uExpected(expectedFES);
    PETSc::Variational::TestFunction  vExpected(expectedFES);
    Problem expected(uExpected, vExpected);
    expected = Integral(gamma * Grad(uExpected), Grad(vExpected))
             - Integral(RealFunction(2.0), vExpected);
    expected.assemble();

    expectSameOwnedVector(
        expected.getLinearSystem().getVector(),
        p.getLinearSystem().getVector());
  }

  // A problem with different global dimensions must produce a different nonzero
  // pattern. Each distinct-size problem uses its own distributed mesh and
  // matrix (the production path for a different mesh); the two patterns must
  // differ in both global dimension and total stored nonzero count.
  TEST(PETSc_MPI_TargetedAssembly, ReassemblyChangesNonzeroPattern)
  {
    const auto& world = *g_world;

    Context::MPI ctx(*g_env, world);

    auto coarseMesh = distributeFromRoot(ctx, 4);
    P1<Real, Mesh<Context::MPI>> coarseFES(coarseMesh);
    PETSc::Variational::TrialFunction uC(coarseFES);
    PETSc::Variational::TestFunction  vC(coarseFES);
    RealFunction gammaC(1.0);
    Problem coarse(uC, vC);
    coarse = Integral(gammaC * Grad(uC), Grad(vC))
           - Integral(RealFunction(1.0), vC);
    coarse.assemble();
    Mat coarseA = coarse.getLinearSystem().getOperator();
    PetscInt rows1 = 0;
    PetscInt cols1 = 0;
    ASSERT_EQ(MatGetSize(coarseA, &rows1, &cols1), PETSC_SUCCESS);
    const PetscInt nz1 = matrixGlobalNonzeros(coarseA);

    auto fineMesh = distributeFromRoot(ctx, 6);
    P1<Real, Mesh<Context::MPI>> fineFES(fineMesh);
    PETSc::Variational::TrialFunction uF(fineFES);
    PETSc::Variational::TestFunction  vF(fineFES);
    RealFunction gammaF(1.0);
    Problem fine(uF, vF);
    fine = Integral(gammaF * Grad(uF), Grad(vF))
         - Integral(RealFunction(1.0), vF);
    fine.assemble();
    Mat fineA = fine.getLinearSystem().getOperator();
    PetscInt rows2 = 0;
    PetscInt cols2 = 0;
    ASSERT_EQ(MatGetSize(fineA, &rows2, &cols2), PETSC_SUCCESS);

    EXPECT_NE(rows1, rows2);
    EXPECT_NE(cols1, cols2);
    EXPECT_NE(nz1, matrixGlobalNonzeros(fineA));
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  [[maybe_unused]] PetscErrorCode ierr =
    PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
