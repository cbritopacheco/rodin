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
  Geometry::MPISharder sharder(mpi);

  if (world.rank() == 0)
  {
    std::cout << "Sharding\n";
    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 32, 32 });
    mesh.getConnectivity().compute(2, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    sharder.shard(partitioner);
    std::cout << "Scatter\n";
    sharder.scatter(0);
  }

  std::cout << "Gather\n";
  auto mesh = sharder.gather(0);

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

  CG(poisson).solve();

  // mpiMesh.getShard().save(
  //     "Gathered" + std::to_string(world.rank()) + ".mesh",
  //     IO::FileFormat::MEDIT);
  // mpiMesh.save("Gathered" + std::to_string(world.rank()) + ".mesh", IO::FileFormat::MEDIT);
}

