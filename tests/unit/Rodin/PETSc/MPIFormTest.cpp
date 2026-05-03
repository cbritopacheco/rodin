/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <functional>

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
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
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

  TEST(PETSc_MPI_Form, DistributedLinearFormUsesMeshCommunicatorAndAssembles)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
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
    EXPECT_EQ(commSize, world.size());

    size_t begin = 0;
    size_t end = 0;
    fes.getOwnershipRange(begin, end);

    PetscInt localSize = 0;
    PetscInt globalSize = 0;
    ierr = VecGetLocalSize(b, &localSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetSize(b, &globalSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(localSize, static_cast<PetscInt>(end - begin));
    EXPECT_EQ(globalSize, static_cast<PetscInt>(fes.getSize()));
  }

  TEST(PETSc_MPI_Form, DistributedBilinearFormUsesMeshCommunicatorAndAssembles)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
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
    EXPECT_EQ(commSize, world.size());

    size_t begin = 0;
    size_t end = 0;
    fes.getOwnershipRange(begin, end);

    PetscInt localRows = 0;
    PetscInt localCols = 0;
    PetscInt globalRows = 0;
    PetscInt globalCols = 0;
    ierr = MatGetLocalSize(A, &localRows, &localCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatGetSize(A, &globalRows, &globalCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    const auto ownedSize = static_cast<PetscInt>(end - begin);
    EXPECT_EQ(localRows, ownedSize);
    EXPECT_EQ(localCols, ownedSize);
    EXPECT_EQ(globalRows, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalCols, static_cast<PetscInt>(fes.getSize()));
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
