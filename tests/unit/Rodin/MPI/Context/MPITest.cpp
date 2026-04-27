#include <gtest/gtest.h>

#include <Rodin/MPI.h>

using namespace Rodin;

static boost::mpi::environment* g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_MPI_Context, ConstGetters)
  {
    const auto& env = *g_env;
    const auto& world = *g_world;

    const Rodin::Context::MPI context(env, world);

    const boost::mpi::communicator& commFromClass = context.getCommunicator();
    const boost::mpi::environment& envFromClass = context.getEnvironment();

    EXPECT_EQ(commFromClass.rank(), world.rank());
    EXPECT_EQ(commFromClass.size(), world.size());
    EXPECT_EQ(&envFromClass, &env);
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
