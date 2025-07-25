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
    std::cout << "Sharding\n";
    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());

    for (auto it = mesh.getCell(); it; ++it)
    {
      mesh.setAttribute(
          {it->getDimension(), it->getIndex()}, partitioner.getPartition(it->getIndex()) + 1);
    }
    mesh.save("Global.mesh", IO::FileFormat::MEDIT);
    sharder.shard(partitioner).scatter(0);
  }

  auto mesh = sharder.gather(0);
  auto& shard = mesh.getShard();
  shard.save("Gathered" + std::to_string(world.rank()) + ".mesh", IO::FileFormat::MEDIT);
  std::cout << "Gathered.\n";
  P1 vh(mesh);
  std::cout << "Created P1 space with " << vh.getSize() << " DOFs.\n";

  PETSc::GridFunction gf(vh);

  // Define problem
  // PETSc::TrialFunction u(vh);
  // PETSc::TestFunction  v(vh);

  // ScalarFunction f = 1;

  // PETSc::Problem poisson(u, v);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(f, v)
  //         + DirichletBC(u, Zero());
  // poisson.assemble();

  // LinearForm lf(v, b);
  // lf = Integral(f, v);
  // lf.assemble();

  // VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // BilinearForm bf(u, v, a);
  // bf = Integral(Grad(u), Grad(v));
  // bf.assemble();

  // poisson.assemble();

  // std::cout << "Matrix :" << "\n";
  // MatView(a, PETSC_VIEWER_STDOUT_WORLD);

  // std::cout << "RHS size: " << "\n";
  // VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // CG(poisson).solve();

  // std::cout << "x after solve:\n";
  // VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  // std::ostringstream mesh_name, sol_name;
  // mesh_name << "mesh." << std::setfill('0') << std::setw(6) << world.rank();
  // sol_name << "sol." << std::setfill('0') << std::setw(6) << world.rank();

  // mesh.save(mesh_name.str(), IO::FileFormat::MFEM);
  // u.getSolution().save(sol_name.str(), IO::FileFormat::MFEM);

  // // mpiMesh.getShard().save(
  // //     "Gathered" + std::to_string(world.rank()) + ".mesh",
  // //     IO::FileFormat::MEDIT);
  // // mpiMesh.save("Gathered" + std::to_string(world.rank()) + ".mesh", IO::FileFormat::MEDIT);

  // PetscFinalize();
}

