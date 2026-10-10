/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include "GaussLobatto.h"

#include <array>
#include <memory>
#include <mutex>
#include <utility>

namespace Rodin::QF
{
  const GaussLobatto& GaussLobatto::get(Geometry::Polytope::Type g, size_t count)
  {
    assert(count >= 2);
    // Match the generic dispatcher's bounded per-thread working set.
    static constexpr size_t HotCapacity = 8;
    static thread_local std::array<const GaussLobatto*, HotCapacity> hot{};
    static thread_local size_t next = 0;
    static FlatMap<std::pair<Geometry::Polytope::Type, size_t>,
      std::unique_ptr<const GaussLobatto>> pool;
    static std::mutex mutex;

    for (const auto* rule : hot)
    {
      if (rule && rule->m_geometry == g && rule->m_nx == count)
        return *rule;
    }

    const GaussLobatto* rule = nullptr;
    {
      const std::lock_guard<std::mutex> lock(mutex);
      const auto key = std::make_pair(g, count);
      auto found = pool.find(key);
      if (found == pool.end())
        found = pool.emplace(key, std::make_unique<GaussLobatto>(g, count)).first;
      rule = found->second.get();
    }
    hot[next] = rule;
    next = (next + 1) % HotCapacity;
    return *rule;
  }
}
