#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <iostream>

#include <Rodin/Context/MPI.h>

namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  Rodin::Context::MPI context(env, world);

  std::cout << "I am process " << context.getWorld().rank() << " of " << context.getWorld().size()
            << "." << std::endl;
  return 0;
}

