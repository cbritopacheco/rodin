/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>

#include <boost/mpi.hpp>

#include <Rodin/Mesh.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  // boost::mpi::environment{argc, argv};
  MPI_Init(&argc, &argv);
  boost::mpi::communicator world;

  Mesh serMesh;
  serMesh.load("../resources/mfem/simple-cantilever2d-example.mesh");

  auto p = serMesh.parallelize(world);
  std::cout << "woof: " << world.rank() << ", " << world.size() << std::endl;
  std::cout << "meow: " << p.getMPIComm().rank() << ", " << p.getMPIComm().size() << std::endl;

  H1 fes(p);

  MPI_Finalize();

  return 0;
}

