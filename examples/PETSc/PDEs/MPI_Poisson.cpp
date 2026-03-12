#include "Rodin/Alert/Info.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Alert/Success.h"
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

static constexpr int ROOT_RANK = 0;

static constexpr int n = 128;


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
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Tetrahedron, { n, n, n });
    Alert::Info() << "Number of cells in the mesh: " << mesh.getCellCount() << Alert::Raise;
    Alert::Info() << "Number of vertices in the mesh: " << mesh.getVertexCount() << Alert::Raise;
    mesh.scale(1.0 / (n - 1));
    Alert::Info() << "Computing mesh connectivity..." << Alert::Raise;
    auto t0 = std::chrono::high_resolution_clock::now();
    const size_t cellDim = mesh.getDimension();
    mesh.getConnectivity().compute(cellDim, cellDim);
    auto t1 = std::chrono::high_resolution_clock::now();

    Alert::Success() << "Computed mesh connectivity in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms." << Alert::Raise;

    Alert::Info() << "Partitioning mesh..." << Alert::Raise;
    t0 = std::chrono::high_resolution_clock::now();
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    t1 = std::chrono::high_resolution_clock::now();

    Alert::Success() << "Partitioned mesh in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms." << Alert::Raise;

    Alert::Info() << "Sharding mesh..." << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    sharder.shard(partitioner);
    t1 = std::chrono::high_resolution_clock::now();

    Alert::Success() << "Sharded mesh in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms." << Alert::Raise;

    Alert::Info() << "Scattering mesh to all processes..." << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    sharder.scatter(ROOT_RANK);
    t1 = std::chrono::high_resolution_clock::now();

    Alert::Success() << "Scattered mesh to all processes in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms." << Alert::Raise;
    Alert::Info() << "Mesh partitioned and scattered to all processes." << Alert::Raise;
  }

  auto mesh = sharder.gather(ROOT_RANK);

  char filename[32];
  std::snprintf(filename, sizeof(filename), "mesh.%06d", world.rank());
  if (world.rank() == ROOT_RANK)
    Alert::Info() << "Saving mesh to " << "mesh.xxxxxx" << "." << Alert::Raise;
  mesh.save(filename);

  if (world.rank() == ROOT_RANK)
    Alert::Info() << "Computing mesh connectivity..." << Alert::Raise;
  mesh.getConnectivity().compute(2, 3);

  if (world.rank() == ROOT_RANK)
    Alert::Info() << "Reconciling." << Alert::Raise;
  mesh.reconcile(2);

  if (world.rank() == ROOT_RANK)
    Alert::Success() << "Saved mesh to " << "mesh.xxxxxx" << "." << Alert::Raise;

  RealFunction f = 1;

  if (world.rank() == ROOT_RANK)
    Alert::Info() << "Constructing finite element space..." << Alert::Raise;

  P1 vh(mesh);

  if (world.rank() == ROOT_RANK)
    Alert::Success() << "Constructed finite element space." << Alert::Raise;

  {
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    // Define problem
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    if (world.rank() == ROOT_RANK)
      Alert::Info() << "Assembling linear system..." << Alert::Raise;

    auto t0 = std::chrono::high_resolution_clock::now();
    poisson.assemble();
    auto t1 = std::chrono::high_resolution_clock::now();

    if (world.rank() == ROOT_RANK)
      Alert::Success() << "Assembled linear system in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                << " ms." << Alert::Raise;

    if (world.rank() == ROOT_RANK)
      Alert::Info() << "Solving linear system..." << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    CG(poisson).solve();
    t1 = std::chrono::high_resolution_clock::now();

    if (world.rank() == ROOT_RANK)
      Alert::Success() << "Solved linear system in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                << " ms." << Alert::Raise;

    std::snprintf(filename, sizeof(filename), "sol.%06d", world.rank());
    u.getSolution().save(filename);
  }

  (void) ierr;
  PetscFinalize();
}

