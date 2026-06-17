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
  Mesh<Context::Local> makeShardableMesh()
  {
    auto mesh =
      Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);

    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh();
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }

    return sharder.gather(0);
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
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

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
