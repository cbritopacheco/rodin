#include <gtest/gtest.h>

#include <Rodin/Context/MPI.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Context_MPI, ConstGetters)
  {
    int argc = 0;
    char** argv = nullptr;
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    const Rodin::Context::MPI context(env, world);

    const boost::mpi::communicator& commFromClass = context.getWorld();
    const boost::mpi::environment& envFromClass = context.getEnvironment();

    EXPECT_EQ(commFromClass.rank(), world.rank());
    EXPECT_EQ(commFromClass.size(), world.size());
    EXPECT_EQ(&envFromClass, &env);
  }
}
