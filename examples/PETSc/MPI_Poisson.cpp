#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>

#include <Rodin/Geometry/BalancedCompactPartitioner.h>

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

  Geometry::MPISharder sharder(mpi);
  if (world.rank() == 0)
  {
    std::cout << "Sharding\n";
    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    sharder.shard(partitioner);
    std::cout << "Scatter\n";
    sharder.scatter(0);
  }

  std::ostringstream mesh_name, sol_name;
  mesh_name << "mesh." << std::setfill('0') << std::setw(6) << world.rank();
  std::string filename = mesh_name.str();

  std::cout << "Gather\n";
  auto mesh = sharder.gather(0);
  std::cout << "Vertex count: " << mesh.getVertexCount() << "\n";

  std::cout << "Local count: " << mesh.getShard().getVertexCount() << "\n";
  std::cout << "Ghosts: " << mesh.getShard().getGhosts()[0].size() << "\n";
  mesh.save(filename, IO::FileFormat::MFEM);
  std::cout << "Saved\n";

  Mat a;
  MatCreate(mpi.getCommunicator(), &a);

  Vec x = PETSC_NULLPTR;

  Vec b;
  VecCreate(mpi.getCommunicator(), &b);

  LinearSystem axb(a, x, b);

  ScalarFunction f = 1;

  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  Problem poisson(u, v, axb);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, Zero());
  poisson.assemble();

  // LinearForm lf(v, b);
  // lf = Integral(f, v);
  // lf.assemble();

  // VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // BilinearForm bf(u, v, a);
  // bf = Integral(Grad(u), Grad(v));
  // bf.assemble();

  // poisson.assemble();

  // CG(poisson).solve();

  // mpiMesh.getShard().save(
  //     "Gathered" + std::to_string(world.rank()) + ".mesh",
  //     IO::FileFormat::MEDIT);
  // mpiMesh.save("Gathered" + std::to_string(world.rank()) + ".mesh", IO::FileFormat::MEDIT);

  PetscFinalize();
}

