#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <iostream>

#include <Rodin/MPI.h>

namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  world.rank();

  Rodin::Context::MPI context(env, world);

  std::cout << "I am process " << context.getCommunicator().rank() << " of " << context.getCommunicator().size()
            << "." << std::endl;
  return 0;
}

