#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI/MPIMesh.h>

namespace mpi = boost::mpi;

int main(int argc, char** argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  std::cout << "I am process " << world.rank() << " of " << world.size()
            << "." << std::endl;

  Rodin::Context::MPI mpi(env, world);
  Rodin::Geometry::MPIMesh mesh(mpi);
}
