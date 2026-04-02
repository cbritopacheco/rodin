/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Implementation of the PolytopeQuadratureFormula dispatcher class.
 */

#include "PolytopeQuadratureFormula.h"

#include <array>
#include <mutex>
#include <shared_mutex>

#include "Centroid.h"
#include "GaussLegendre.h"
#include "GrundmannMoller.h"

namespace Rodin::QF
{
  namespace
  {
    struct Key
    {
      size_t order;
      Geometry::Polytope::Type g;

      friend bool operator<(const Key& a, const Key& b)
      {
        if (a.g != b.g)
          return a.g < b.g;
        return a.order < b.order;
      }

      friend bool operator==(const Key& a, const Key& b)
      {
        return a.order == b.order && a.g == b.g;
      }
    };

    /**
     * @brief Tiny per-thread hot cache for canonical quadrature pointers.
     *
     * This cache does not own anything. It only mirrors pointers owned by the
     * global canonical pool below. Entries are replaced in round-robin order.
     */
    struct HotCache
    {
      static constexpr size_t Capacity = 8;

      struct Entry
      {
        Key key{0, Geometry::Polytope::Type::Point};
        const QuadratureFormulaBase* ptr = nullptr;
        bool valid = false;
      };

      std::array<Entry, Capacity> entries;
      size_t next = 0;
    };

    /**
     * @brief Process-wide canonical pool of generic quadrature formulas.
     *
     * This pool owns one unique quadrature formula object per `(order, geometry)`
     * pair. Returned references remain valid until program termination.
     */
    struct CanonicalPool
    {
      FlatMap<Key, std::unique_ptr<const QuadratureFormulaBase>> formulas;
      std::mutex mutex;
    };
  }

  std::unique_ptr<const QuadratureFormulaBase>
  PolytopeQuadratureFormula::build(size_t order, Geometry::Polytope::Type g)
  {
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
        return std::make_unique<Centroid>(g);

      case Geometry::Polytope::Type::Segment:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n);
      }

      case Geometry::Polytope::Type::Triangle:
      {
        const size_t s = (order + 1) / 2;
        return std::make_unique<GrundmannMoller>(s, g);
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        const size_t i = (order / 2) * 2 + 1;
        return std::make_unique<GrundmannMoller>(i / 2, g);
      }

      case Geometry::Polytope::Type::Wedge:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n, n);
      }

      default:
      {
        assert(false);
        return {};
      }
    }
  }

  const QuadratureFormulaBase&
  PolytopeQuadratureFormula::get(size_t order, Geometry::Polytope::Type g)
  {
    static CanonicalPool pool;
    static thread_local HotCache hot;

    const Key key{order, g};

    // ------------------------------------------------------------------------
    // Tier 1: tiny thread-local hot cache
    // ------------------------------------------------------------------------
    for (auto& e : hot.entries)
    {
      if (e.valid && e.key == key)
      {
        assert(e.ptr);
        return *e.ptr;
      }
    }

    // ------------------------------------------------------------------------
    // Tier 2: process-wide canonical pool
    // ------------------------------------------------------------------------
    const QuadratureFormulaBase* ptr = nullptr;

    {
      std::lock_guard<std::mutex> lock(pool.mutex);

      auto it = pool.formulas.find(key);
      if (it == pool.formulas.end())
        it = pool.formulas.emplace(key, build(order, g)).first;

      assert(it->second);
      ptr = it->second.get();
    }

    assert(ptr);

    // ------------------------------------------------------------------------
    // Update thread-local hot cache (round-robin replacement)
    // ------------------------------------------------------------------------
    auto& e = hot.entries[hot.next];
    e.key = key;
    e.ptr = ptr;
    e.valid = true;
    hot.next = (hot.next + 1) % HotCache::Capacity;

    return *ptr;
  }

  PolytopeQuadratureFormula::PolytopeQuadratureFormula(
      size_t order,
      Geometry::Polytope::Type g)
    : m_qf(build(order, g)),
      m_order(order),
      m_geometry(g)
  {}
}
