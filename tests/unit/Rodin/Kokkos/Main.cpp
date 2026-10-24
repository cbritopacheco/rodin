/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

/**
 * @brief Runs tests within the lifetime of the Kokkos runtime.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @returns GoogleTest exit status.
 */
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  ::Kokkos::ScopeGuard runtime(argc, argv);
  return RUN_ALL_TESTS();
}
