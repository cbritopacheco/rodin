/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_REVISION_H
#define RODIN_REVISION_H

/**
 * @file
 * @brief Defines a revision stamp for cache validation.
 */

#include <atomic>
#include <cstdint>

namespace Rodin
{
  /**
   * @brief Copyable snapshot of an object's revision.
   *
   * Issued stamps are distinct within the process, including across owners.
   * The owner stores a fresh stamp when the state tracked by a cache changes.
   */
  class Revision
  {
    public:
      /// @brief Constructs the initial, unissued stamp.
      constexpr Revision() noexcept
        : m_value(0)
      {}

      /// @returns A fresh stamp, leaving this snapshot intact.
      Revision next() const noexcept
      {
        static std::atomic<std::uint64_t> counter{0};
        return Revision(counter.fetch_add(1, std::memory_order_relaxed) + 1);
      }

      /// @param other Stamp to compare. @returns Whether both stamps identify the same revision.
      constexpr bool operator==(const Revision& other) const noexcept = default;

      /// @returns Numeric stamp for interoperability with integer revision keys.
      explicit constexpr operator std::uint64_t() const noexcept
      {
        return m_value;
      }

    private:
      /// @param value Fresh stamp issued by next().
      explicit constexpr Revision(std::uint64_t value) noexcept
        : m_value(value)
      {}

      std::uint64_t m_value;
  };
}

#endif
