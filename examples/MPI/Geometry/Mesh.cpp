#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

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
  Geometry::MPISharder sharder(mpi);
  if (world.rank() == 0)
  {

    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 2);
    Geometry::BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(world.size());
    for (auto it = mesh.getCell(); it; ++it)
    {
      mesh.setAttribute({it->getDimension(), it->getIndex()}, partitioner.getPartition(it->getIndex()) + 1);
    }
    mesh.save("Global.mesh", IO::FileFormat::MEDIT);
    std::cout << "Sharding into " << partitioner.getCount() << "...\n";

    sharder.shard(partitioner);

    std::cout << "CellCount: " << sharder.getShard(0).getCellCount() << std::endl;
    std::cout << "FaceCount: " << sharder.getShard(0).getFaceCount() << std::endl;
    std::cout << "VertexCount: " << sharder.getShard(0).getVertexCount() << std::endl;
    std::cout << "Flags: " << sharder.getShard(0).getFlags(0).size() << std::endl;

    std::cout << "Shard 1\n";
    for (auto it = sharder.getShard(0).getCell(); it; ++it)
    {
      if (sharder.getShard(0).isGhost(it->getDimension(), it->getIndex()))
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 4);
      else if (sharder.getShard(0).isOwned(it->getDimension(), it->getIndex()))
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 8);
    }

    for (auto it = sharder.getShard(0).getFace(); it; ++it)
    {
      if (sharder.getShard(0).isGhost(it->getDimension(), it->getIndex()))
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 5);
      else if (sharder.getShard(0).isOwned(it->getDimension(), it->getIndex()))
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 10);
    }

    for (auto it = sharder.getShard(0).getVertex(); it; ++it)
    {
      if (sharder.getShard(0).isGhost(it->getDimension(), it->getIndex()))
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 6);
      else if (sharder.getShard(0).isOwned(it->getDimension(), it->getIndex()))
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 12);
    }

    sharder.getShard(0).save("Shard.0.mesh", IO::FileFormat::MEDIT);

    std::cout << "Shard 2\n";
    for (auto it = sharder.getShard(1).getCell(); it; ++it)
    {
      if (sharder.getShard(1).isGhost(it->getDimension(), it->getIndex()))
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 7);
      else if (sharder.getShard(1).isOwned(it->getDimension(), it->getIndex()))
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 14);
    }

    for (auto it = sharder.getShard(1).getFace(); it; ++it)
    {
      if (sharder.getShard(1).isGhost(it->getDimension(), it->getIndex()))
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 8);
      else if (sharder.getShard(1).isOwned(it->getDimension(), it->getIndex()))
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 16);
    }

    for (auto it = sharder.getShard(1).getVertex(); it; ++it)
    {
      if (sharder.getShard(1).isGhost(it->getDimension(), it->getIndex()))
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 9);
      else if (sharder.getShard(1).isOwned(it->getDimension(), it->getIndex()))
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 18);
    }
    sharder.getShard(1).save("Shard.1.mesh", IO::FileFormat::MEDIT);

    std::cout << "Scattering\n";
    sharder.scatter(0);
  }

  std::cout << "Gathering\n";
  auto mpiMesh = sharder.gather(0);
  mpiMesh.getShard().save(
      "Gathered" + std::to_string(world.rank()) + ".mesh",
      IO::FileFormat::MEDIT);
  // mpiMesh.save("Gathered" + std::to_string(world.rank()) + ".mesh", IO::FileFormat::MEDIT);
}
