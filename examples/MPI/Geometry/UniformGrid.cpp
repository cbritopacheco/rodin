/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/mpi/environment.hpp>

#include <Rodin/Alert.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Geometry.h>

#include <Rodin/MPI.h>

using namespace Rodin;
using namespace Geometry;

namespace mpi = boost::mpi;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  constexpr size_t n = 128;
  constexpr Geometry::Polytope::Type g = Geometry::Polytope::Type::Tetrahedron;

  Context::MPI mpi(env, world);
  auto t0 = std::chrono::high_resolution_clock::now();
  auto mesh = Mesh<Context::MPI>::UniformGrid(mpi, g, { n, n, n });
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = t1 - t0;

  size_t nv = mesh.getVertexCount();
  size_t nc = mesh.getCellCount();

  if (world.rank() == 0)
  {
    Alert::Info()
      << "Generated uniform grid with " << nv << " vertices and "
      << nc << " cells in " << elapsed.count() << " seconds."
      << Alert::Raise;
  }

  IO::XDMF xdmf(world, "UniformGrid");
  xdmf.grid().setMesh(mesh);
  xdmf.write();
  xdmf.close();

  return 0;
}




