/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <type_traits>

#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  static_assert(std::is_base_of_v<Context::Local, TestContext>);
  static_assert(std::is_same_v<TestContext::DistributionType, Context::Local>);
  static_assert(std::is_same_v<TestContext::BackendType, TestBackend>);
  static_assert(std::is_same_v<TestMemorySpace, TestExecutionSpace::memory_space>);
#ifdef KOKKOS_ENABLE_CUDA
  static_assert(std::is_same_v<TestExecutionSpace, ::Kokkos::Cuda>);
#endif

  /// @brief A copied context retains a usable local execution association.
  void checkCopyRetainsExecution()
  {
    TestContext original{Context::Local{}, TestExecutionSpace{}};
    TestContext copied(original);
    EXPECT_EQ(&original.getDistribution(), static_cast<const Context::Local*>(&original));
    EXPECT_EQ(&copied.getDistribution(), static_cast<const Context::Local*>(&copied));
    ::Kokkos::View<Index*, TestMemorySpace> values("context test values", 17);
    ::Kokkos::parallel_for(
      "context test execution",
      ::Kokkos::RangePolicy<TestExecutionSpace>(copied.getExecutionSpace(), 0, 17),
      KOKKOS_LAMBDA(Index i) { values(i) = i * i; });
    auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, values);
    for (Index i = 0; i < 17; ++i)
      EXPECT_EQ(host(i), i * i);
  }

  /// @brief Exercises CopyRetainsExecution with CUDA-compatible kernel scope.
  TEST(Kokkos_Context, CopyRetainsExecution)
  {
    checkCopyRetainsExecution();
  }
}
