#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI/MPIMesh.h>
#include <Rodin/MPI/MPISharder.h>

#include <Rodin/Geometry/BalancedCompactPartitioner.h>

namespace mpi = boost::mpi;

using namespace Rodin;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  std::cout << "I am process " << world.rank() << " of " << world.size()
            << "." << std::endl;

  Context::MPI mpi(env, world);
  Geometry::MPISharder sharder(mpi);
  if (world.rank() == 0)
  {

    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 32, 32 });
    mesh.getConnectivity().compute(2, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    for (auto it = mesh.getCell(); it; ++it)
    {
      mesh.setAttribute({it->getDimension(), it->getIndex()}, partitioner.getPartition(it->getIndex()) + 1);
    }
    mesh.save("Global.mesh");
    std::cout << "Sharding into " << partitioner.getCount() << "...\n";
    sharder.shard(partitioner);
    std::cout << "miaow: " << sharder.getShard(0).getCellCount() << std::endl;
    sharder.getShard(0).save("Shard.0.mesh");
    sharder.getShard(1).save("Shard.1.mesh");
    std::cout << "Scattering\n";
    // sharder.scatter(0);
  }

  std::cout << "Gathering\n";
  // auto mpiMesh = sharder.gather(0);
}
