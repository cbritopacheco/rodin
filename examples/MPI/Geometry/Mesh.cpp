#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>

#include <Rodin/Geometry/BalancedCompactPartitioner.h>

namespace mpi = boost::mpi;

using namespace Rodin;
using namespace Rodin::MPI;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  Context::MPI mpi(env, world);
  Sharder sharder(mpi);
  if (world.rank() == 0)
  {

    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 32, 32 });
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(2, 1);
    std::cout << mesh.getFaceCount() << std::endl;
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    for (auto it = mesh.getCell(); it; ++it)
    {
      mesh.setAttribute(
          {it->getDimension(), it->getIndex()}, partitioner.getPartition(it->getIndex()) * 3 + 1);
    }
    mesh.save("Global.mesh", IO::FileFormat::MEDIT);
    std::cout << "Sharding into " << partitioner.getCount() << "...\n";

    sharder.shard(partitioner);
    auto& shards = sharder.getShards();

    std::cout << "Shard 1\n";
    for (size_t s = 0; s < world.size(); s++)
    {
      for (auto it = shards[s].getCell(); it; ++it)
      {
        if (shards[s].isGhost(it->getDimension(), it->getIndex()))
          shards[s].setAttribute({ it->getDimension(), it->getIndex() }, 4);
        else if (shards[s].isOwned(it->getDimension(), it->getIndex()))
          shards[s].setAttribute({ it->getDimension(), it->getIndex() }, 8);
      }

      for (auto it = shards[s].getFace(); it; ++it)
      {
        if (shards[s].isGhost(it->getDimension(), it->getIndex()))
          shards[s].setAttribute({ it->getDimension(), it->getIndex() }, 4);
        else if (shards[s].isOwned(it->getDimension(), it->getIndex()))
          shards[s].setAttribute({ it->getDimension(), it->getIndex() }, 8);
      }

      for (auto it = shards[s].getVertex(); it; ++it)
      {
        if (shards[s].isGhost(it->getDimension(), it->getIndex()))
          shards[s].setAttribute({ it->getDimension(), it->getIndex() }, 4);
        else if (shards[s].isOwned(it->getDimension(), it->getIndex()))
          shards[s].setAttribute({ it->getDimension(), it->getIndex() }, 8);
      }
      shards[s].save("Shard." + std::to_string(s) + ".mesh", IO::FileFormat::MEDIT);
    }

    std::cout << "Scattering\n";
    sharder.scatter(0);
  }

  std::cout << "Gathering\n";
  auto mpiMesh = sharder.gather(0);
  mpiMesh.getShard().save(
      "Gathered." + std::to_string(world.rank()) + ".mesh",
      IO::FileFormat::MEDIT);
}
