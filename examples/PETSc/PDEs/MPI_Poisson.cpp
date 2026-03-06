#include "Rodin/Alert/Raise.h"
#include "Rodin/PETSc/Variational/GridFunction.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Alert.h>

namespace mpi = boost::mpi;

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

static constexpr Index ROOT_RANK = 0;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  Context::MPI mpi(env, world);

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  assert(ierr == PETSC_SUCCESS);

  Rodin::MPI::Sharder sharder(mpi);
  if (world.rank() == ROOT_RANK)
  {
    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Hexahedron, { 64, 64, 64 });
    mesh.scale(1.0 / 63.0);
    Alert::Info() << "Computing mesh connectivity..." << Alert::Raise;
    const size_t cellDim = mesh.getDimension();
    mesh.getConnectivity().compute(cellDim, cellDim);
    mesh.getConnectivity().compute(cellDim - 1, cellDim);
    Alert::Info() << "Partitioning mesh..." << Alert::Raise;
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    Alert::Info() << "Sharding mesh..." << Alert::Raise;
    sharder.shard(partitioner);
    Alert::Info() << "Scattering mesh to all processes..." << Alert::Raise;
    sharder.scatter(ROOT_RANK);
    Alert::Info() << "Mesh partitioned and scattered to all processes." << Alert::Raise;
  }

  PetscPrintf(world, "Gathering mesh on root process for visualization.\n");
  auto mesh = sharder.gather(ROOT_RANK);
  PetscPrintf(world, "Gathered mesh on root process.\n");

  char filename[32];
  std::snprintf(filename, sizeof(filename), "mesh.%06d", world.rank());
  PetscPrintf(world, "Saving mesh to %s.\n", filename);
  mesh.save(filename);
  PetscPrintf(world, "Saved mesh to %s.\n", filename);

  RealFunction f = 1;

  PetscPrintf(world, "Constructing finite element space.\n");
  P1 vh(mesh);
  PetscPrintf(world, "Constructed finite element space.\n");

  {
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    // Define problem
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    PetscPrintf(world, "Assembling.\n");
    poisson.assemble();
    PetscPrintf(world, "Assembled linear system.\n");
    PetscPrintf(world, "Solving.\n");
    CG(poisson).solve();
    PetscPrintf(world, "Solved linear system.\n");

    std::snprintf(filename, sizeof(filename), "sol.%06d", world.rank());
    u.getSolution().save(filename);
  }

  PetscFinalize();
}

