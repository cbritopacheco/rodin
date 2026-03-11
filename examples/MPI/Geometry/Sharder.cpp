#include "Rodin/Alert/Success.h"
#include "Rodin/IO/ForwardDecls.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/Alert/Info.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>

#include <Rodin/Geometry/BalancedCompactPartitioner.h>

namespace mpi = boost::mpi;

using namespace Rodin;

static constexpr int ROOT_RANK = 0;

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
    mesh = mesh.UniformGrid(Geometry::Polytope::Type::Tetrahedron, { 32, 32, 32 });
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    Alert::Info() << "Mesh generation time: " << duration << " ms" << Alert::Raise;

    Alert::Info() << "Number of cells in the mesh: " << mesh.getCellCount() << Alert::Raise;
    Alert::Info() << "Number of vertices in the mesh: " << mesh.getVertexCount() << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    Alert::Info() << "Computing mesh connectivity..." << Alert::Raise;
    const size_t cellDim = mesh.getDimension();
    mesh.getConnectivity().compute(cellDim, cellDim);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    Alert::Info() << "Computed mesh connectivity in " << duration << " ms." << Alert::Raise;

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
    sharder.scatter(ROOT_RANK);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    Alert::Info() << "Scattering time: " << duration << " ms" << Alert::Raise;
  }

  if (world.rank() == ROOT_RANK)
    Alert::Info() << "Gathering mesh on all processes..." << Alert::Raise;

  auto mesh = sharder.gather(0);

  mesh.getShard().save("Shard." + std::to_string(world.rank()) + ".mesh", IO::FileFormat::MEDIT);

  if (world.rank() == ROOT_RANK)
    Alert::Success() << "Saved shard meshes as Shard.x.mesh (MEDIT)" << Alert::Raise;
}

