/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "HilbertOrdering.h"

namespace Rodin::Geometry
{
  void HilbertOrdering::transform(std::vector<std::uint64_t>& coordinates, size_t bits)
  {
    // Exchange/invert coordinate blocks to undo the Hilbert orientation.
    const std::uint64_t highest = std::uint64_t(1) << (bits - 1);
    for (std::uint64_t q = highest; q > 1; q >>= 1)
    {
      const std::uint64_t mask = q - 1;
      for (size_t i = 0; i < coordinates.size(); ++i)
      {
        if (coordinates[i] & q)
          coordinates[0] ^= mask;
        else
        {
          const std::uint64_t exchange = (coordinates[0] ^ coordinates[i]) & mask;
          coordinates[0] ^= exchange;
          coordinates[i] ^= exchange;
        }
      }
    }
    // Gray encode and apply the final orientation correction.
    for (size_t i = 1; i < coordinates.size(); ++i)
      coordinates[i] ^= coordinates[i - 1];
    std::uint64_t correction = 0;
    for (std::uint64_t q = highest; q > 1; q >>= 1)
    {
      if (coordinates.back() & q)
        correction ^= q - 1;
    }
    for (auto& coordinate : coordinates)
      coordinate ^= correction;
  }
}
