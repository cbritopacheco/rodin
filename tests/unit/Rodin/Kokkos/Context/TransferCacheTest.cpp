/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "Rodin/Kokkos/Utility.h"

#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  using TransferCache = Rodin::Kokkos::Internal::TransferCache<TestMemorySpace>;

  /// @brief Host values distinguished by a shift.
  std::vector<Real> makeValues(size_t size, Real shift)
  {
    std::vector<Real> res(size);
    for (size_t i = 0; i < size; ++i)
      res[i] = shift + Real(0.25) * i;
    return res;
  }

  /// @brief Reads a resident copy back to the host.
  std::vector<Real> readBack(const TransferCache::Array& array)
  {
    const auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, array);
    return std::vector<Real>(host.data(), host.data() + host.extent(0));
  }

  /// @brief Copies are shared by content and never overwritten while referenced.
  void checkTransferCache()
  {
    const TestExecutionSpace execution;
    TransferCache cache;
    std::vector<const Real*> addresses;
    {
      const auto a = cache.upload(execution, makeValues(32, 0), "a");
      EXPECT_EQ(readBack(a), makeValues(32, 0));
      const auto same = cache.upload(execution, makeValues(32, 0), "a");
      EXPECT_EQ(same.data(), a.data());

      // Every kept copy is referenced: the next contents go elsewhere.
      std::vector<TransferCache::Array> held{a};
      for (size_t k = 1; k < TransferCache::Capacity + 2; ++k)
        held.push_back(cache.upload(execution, makeValues(32, Real(k)), "b"));
      for (size_t k = 0; k < held.size(); ++k)
      {
        EXPECT_EQ(readBack(held[k]), makeValues(32, Real(k))) << k;
        for (size_t j = 0; j < k; ++j)
          EXPECT_NE(held[k].data(), held[j].data()) << k << " " << j;
      }
      for (size_t k = 0; k < TransferCache::Capacity; ++k)
        addresses.push_back(held[k].data());
    }

    // Released copies take the next contents in place.
    const auto next = cache.upload(execution, makeValues(32, 100), "c");
    EXPECT_EQ(readBack(next), makeValues(32, 100));
    EXPECT_NE(
      std::find(addresses.begin(), addresses.end(), next.data()), addresses.end());
    const auto kept = cache.upload(execution, makeValues(32, 3), "c");
    EXPECT_EQ(readBack(kept), makeValues(32, 3));
    EXPECT_EQ(readBack(next), makeValues(32, 100));

    // Another size replaces a released copy.
    const auto other = cache.upload(execution, makeValues(7, 5), "d");
    EXPECT_EQ(readBack(other), makeValues(7, 5));
    EXPECT_EQ(readBack(next), makeValues(32, 100));
    EXPECT_EQ(readBack(kept), makeValues(32, 3));
  }

  TEST(Kokkos_TransferCache, SharesByContentAndKeepsReferencedCopies)
  {
    checkTransferCache();
  }
}
