#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/MPI/Variational/P1.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  Context::MPI mpi(env, world);

  Mesh<Context::MPI> mesh(mpi);
  P1 fes(mesh);
}

