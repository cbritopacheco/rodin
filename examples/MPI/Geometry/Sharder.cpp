#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/Alert/Info.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>

#include <Rodin/Geometry/BalancedCompactPartitioner.h>

namespace mpi = boost::mpi;

using namespace Rodin;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  Context::MPI mpi(env, world);
  Rodin::MPI::Sharder sharder(mpi);

  if (world.rank() == 0)
  {
    Geometry::LocalMesh mesh;
    auto t0 = std::chrono::high_resolution_clock::now();
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Tetrahedron, { 64, 64, 64 });
    mesh.getConnectivity().compute(2, 2);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    Alert::Info() << "Mesh generation time: " << duration << " ms" << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    Alert::Info() << "Partitioning time: " << duration << " ms" << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    sharder.shard(partitioner);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    Alert::Info() << "Sharding time: " << duration << " ms" << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    sharder.scatter(0);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    Alert::Info() << "Scattering time: " << duration << " ms" << Alert::Raise;
  }

  Geometry::MPIMesh mesh = sharder.gather(0);

  mesh.getShard().save("Shard" + std::to_string(world.rank()) + ".mesh");
}

