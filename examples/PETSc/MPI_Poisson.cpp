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
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 2);
    mesh.save("Poisson.mesh");
    // mesh.getConnectivity().compute(1, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    sharder.shard(partitioner).scatter(ROOT_RANK);
  }

  auto mesh = sharder.gather(ROOT_RANK);
  P1 vh(mesh);

  mesh.save("Poisson." + std::to_string(world.rank()) + ".mesh");

  ScalarFunction f = 1;

  {
    PETSc::Variational::GridFunction gf(vh);
    gf = Cos(F::x);
    gf.save("Poisson." + std::to_string(world.rank()) + ".gf");
  }

  // Vec v;
  // VecCreate(PETSC_COMM_WORLD, &v);

  // size_t sz;
  // if (world.rank() == 0)
  // {
  //   sz = 7;
  // }
  // else
  // {
  //   sz = 8;
  // }
  // VecSetSizes(v, sz, 10);
  // VecSetFromOptions(v);

  // std::vector<int> ghosts;
  // if (world.rank() == 0)
  // {
  //   ghosts = { 2, 3 };
  // }
  // else
  // {
  //   ghosts = { 4, 5, 6 };
  // }

  // VecMPISetGhost(v, ghosts.size(), ghosts.data());


  PetscFinalize();
}

