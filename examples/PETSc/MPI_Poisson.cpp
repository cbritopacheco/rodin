#include "Rodin/Variational/ForwardDecls.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>

#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <string>

namespace mpi = boost::mpi;

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  Context::MPI mpi(env, world);

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  assert(ierr == PETSC_SUCCESS);

  Rodin::MPI::Sharder sharder(mpi);
  if (world.rank() == 0)
  {
    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    sharder.shard(partitioner).scatter(0);
  }

  auto mesh = sharder.gather(0);
  P1 vh(mesh);

  {
    PETSc::GridFunction gf(vh);
  }

  PetscFinalize();
}

