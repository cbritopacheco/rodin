#include "Rodin/Variational/ForwardDecls.h"
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
  sleep(5);
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
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    sharder.shard(partitioner).scatter(0);
  }

  // std::cout << "Gather\n";
  auto mesh = sharder.gather(0);
  // mesh.getShard().getConnectivity().compute(1, 2);
  // std::cout << "Vertex count: " << mesh.getVertexCount() << "\n";
  // std::cout << "Shard: " << world.rank() <<  " Local count: "
  // << mesh.getShard().getVertexCount() << "\n";

  auto& shard = mesh.getShard();

  mesh.getShard().save(
      "Gathered" + std::to_string(world.rank()) + ".mesh",
      IO::FileFormat::MEDIT);
  // print polytope map
  for (size_t d = 0; d < mesh.getDimension(); ++d)
  {
    std::cout << "Dimension: " << d << "\n";
    const auto& map = shard.getPolytopeMap(d);
    for (const auto& [local, global] : map.left)
    {
      std::cout << "Local: " << local << ", Global: " << global << "\n";
    }
  }

  // Print halo
  if (world.rank() == 0)
  {
    std::cout << "Halo information:\n";
    for (size_t d = 0; d <= mesh.getDimension(); ++d)
    {
      std::cout << "Dimension: " << d << "\n";
      for (Index i = 0; i < shard.getPolytopeCount(d); ++i)
      {
          std::cout << "Owned polytope: " << i << "\n";
          const auto& halo = shard.getHalo(d, i);
          std::cout << "Halo: ";
          for (const auto& h : halo)
          {
            std::cout << h << " ";
          }
          std::cout << "\n";
        }
    }
  }

  // P1 vh(mesh);
  // PETSc::GridFunction gf(vh);

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

