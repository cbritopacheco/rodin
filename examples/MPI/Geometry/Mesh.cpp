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

  std::cout << "I am process " << world.rank() << " of " << world.size()
            << "." << std::endl;

  Context::MPI mpi(env, world);
  Geometry::MPISharder sharder(mpi);
  if (world.rank() == 0)
  {

    Geometry::LocalMesh mesh;
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, { 64, 64 });
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
    std::cout << "miaow: " << sharder.getShard(0).getCellCount() << std::endl;

    for (auto it = sharder.getShard(0).getCell(); it; ++it)
    {
      if (sharder.getShard(0).isGhost(it->getDimension(), it->getIndex()))
      {
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 1);
      }
      else
      {
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 2);
      }
    }

    for (auto it = sharder.getShard(1).getCell(); it; ++it)
    {
      if (sharder.getShard(1).isGhost(it->getDimension(), it->getIndex()))
      {
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 2);
      }
      else
      {
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 1);
      }
    }

    for (auto it = sharder.getShard(0).getVertex(); it; ++it)
    {
      if (sharder.getShard(0).isGhost(it->getDimension(), it->getIndex()))
      {
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 5);
      }
      else
      {
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 10);
      }
    }

    for (auto it = sharder.getShard(1).getVertex(); it; ++it)
    {
      if (sharder.getShard(1).isGhost(it->getDimension(), it->getIndex()))
      {
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 5);
      }
      else
      {
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 10);
      }
    }

    for (auto it = sharder.getShard(0).getFace(); it; ++it)
    {
      sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 2);
      if (sharder.getShard(0).isGhost(it->getDimension(), it->getIndex()))
      {
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 4);
      }
      else
      {
        sharder.getShard(0).setAttribute({ it->getDimension(), it->getIndex() }, 12);
      }
    }

    for (auto it = sharder.getShard(1).getFace(); it; ++it)
    {
      if (sharder.getShard(1).isGhost(it->getDimension(), it->getIndex()))
      {
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 4);
      }
      else
      {
        sharder.getShard(1).setAttribute({ it->getDimension(), it->getIndex() }, 12);
      }
    }

    sharder.getShard(0).save("Shard.0.mesh", IO::FileFormat::MEDIT);
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
