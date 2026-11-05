/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file Main.cpp
 * @brief Entry point of the distributed device tests.
 *
 * The Poisson and constraint sources instantiate the same device assembly
 * templates; linking them in one executable checks that their kernels do not
 * depend on the translation unit that launches them.
 */

#include "Common.h"

/**
 * @brief Runs tests within the lifetimes of the MPI, Kokkos and PETSc runtimes.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @returns GoogleTest exit status.
 */
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  Rodin::Tests::Unit::g_env = &env;
  Rodin::Tests::Unit::g_world = &world;
  ::testing::InitGoogleTest(&argc, argv);
  int res = 0;
  {
    ::Kokkos::ScopeGuard runtime(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    res = RUN_ALL_TESTS();
    PetscFinalize();
  }
  return res;
}
