/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPEQUADRATURE_H
#define RODIN_GEOMETRY_POLYTOPEQUADRATURE_H

/**
 * @file
 * @brief Cached quadrature data attached to mesh polytopes.
 *
 * This file defines:
 * - @ref Rodin::Geometry::PolytopeQuadrature, which stores the mapped
 *   quadrature points associated with a specific quadrature formula on a
 *   specific polytope.
 * - @ref Rodin::Geometry::PolytopeQuadratureIndex, a mesh-owned, thread-safe
 *   cache indexed by polytope dimension, polytope index, and quadrature
 *   formula identity.
 *
 * # Purpose
 *
 * Reference quadrature formulas belong to the @ref Rodin::QF module and live on
 * reference polytopes. In contrast, this file provides a geometry-side cache of
 * their realization on concrete mesh polytopes.
 *
 * For a polytope @f$ \tau @f$ and a reference quadrature formula
 * @f$ \{ (\hat x_q, w_q) \}_{q=1}^{n_q} @f$, the realized quadrature consists
 * of the mapped geometric points:
 *
 * @f[
 *   p_q = x_\tau(\hat x_q), \qquad q = 1, \dots, n_q,
 * @f]
 *
 * where @f$ x_\tau @f$ is the polytope transformation.
 *
 * # Cache organization
 *
 * The cache is organized as:
 *
 * @f[
 *   (d, i) \mapsto \text{Slot}_{d,i}
 * @f]
 *
 * where:
 * - @f$ d @f$ is the polytope dimension,
 * - @f$ i @f$ is the polytope index within that dimension.
 *
 * Each slot stores:
 * - one hot entry for the fastest repeated-hit access path,
 * - one small bounded ring buffer of owned cached quadratures.
 *
 * Lookup order is:
 * - hot entry,
 * - bounded per-slot cache,
 * - lazy construction on miss.
 *
 * Hits in the bounded cache are promoted to the hot slot. Misses are inserted
 * into the bounded cache in round-robin order and then promoted to the hot slot.
 *
 * # Key
 *
 * The cache key is the address of the quadrature formula object:
 *
 * @f[
 *   \texttt{const QF::QuadratureFormulaBase*}
 * @f]
 *
 * This design assumes that callers obtain quadrature formulas from a canonical
 * source such as @ref Rodin::QF::PolytopeQuadratureFormula::get, so that the
 * same logical formula has a stable object identity.
 */

#include <cassert>
#include <utility>
#include <array>
#include <deque>
#include <memory>
#include <mutex>
#include <vector>
#include <shared_mutex>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Types.h"
#include "Rodin/QF/ForwardDecls.h"
#include "Point.h"

namespace Rodin::Geometry
{
  /**
   * @brief Cached quadrature attached to one specific polytope.
   *
   * A PolytopeQuadrature stores the mapped quadrature points corresponding to a
   * single quadrature formula on a single concrete polytope. It is:
   * - not a reference quadrature formula,
   * - not a variational quadrature rule,
   * - not a container of weights.
   *
   * Its purpose is to cache the geometric realization of a quadrature formula
   * on a polytope so that repeated integration passes can reuse mapped
   * @ref Rodin::Geometry::Point objects.
   */
  class PolytopeQuadrature
  {
    public:
      /**
       * @brief Constructs the realized quadrature on a polytope.
       * @param[in] polytope Concrete mesh polytope
       * @param[in] qf Reference quadrature formula
       */
      PolytopeQuadrature(
          const Polytope& polytope,
          const QF::QuadratureFormulaBase& qf);

      /**
       * @brief Gets the reference quadrature formula used to build this object.
       * @returns Reference to the quadrature formula
       */
      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        assert(m_qf);
        return *m_qf;
      }

      /**
       * @brief Gets the number of mapped quadrature points.
       * @returns Number of quadrature points
       */
      size_t getSize() const
      {
        return m_ps.size();
      }

      /**
       * @brief Gets a mapped quadrature point by index.
       * @param[in] i Quadrature point index
       * @returns Reference to the mapped point
       */
      const Point& getPoint(size_t i) const
      {
        return m_ps.at(i);
      }

    private:
      const QF::QuadratureFormulaBase* m_qf = nullptr; ///< Source quadrature formula
      std::vector<Point> m_ps;                         ///< Mapped quadrature points
  };

  /**
   * @brief Thread-safe cache of polytope quadratures.
   *
   * This class maintains a mesh-owned cache of @ref PolytopeQuadrature objects,
   * indexed by:
   * - polytope dimension,
   * - polytope index,
   * - quadrature formula identity.
   *
   * # Slot design
   *
   * Each polytope slot contains:
   * - a hot entry storing a non-owning mirror of the most recently promoted
   *   cached quadrature,
   * - a fixed-capacity array of owned entries used as a bounded ring buffer.
   *
   * The hot entry exists only to accelerate repeated hits. Ownership remains
   * solely in the bounded cache entries.
   *
   * # Eviction policy
   *
   * The bounded cache uses round-robin replacement once full. This is
   * appropriate because:
   * - the per-polytope cache is intentionally small,
   * - the lookup cost is linear in a tiny fixed capacity,
   * - the hot entry captures the most common repeated access pattern.
   *
   * # Thread safety
   *
   * - Each dimension bucket owns a shared mutex.
   * - Each polytope slot owns its own shared mutex.
   * - Reads first take shared locks.
   * - Misses upgrade to exclusive slot access only where insertion is required.
   *
   * # Lifecycle
   *
   * This cache is derived state. It should be cleared whenever the mesh
   * geometry changes.
   *
   * @note This cache is intentionally transient and is not serializable.
   */
  class PolytopeQuadratureIndex
  {
    public:
      /**
       * @brief Default constructor.
       */
      PolytopeQuadratureIndex() = default;

      /**
       * @brief Destructor.
       */
      ~PolytopeQuadratureIndex() = default;

      /**
       * @brief Copy constructor (deleted).
       */
      PolytopeQuadratureIndex(const PolytopeQuadratureIndex&) = delete;

      /**
       * @brief Copy assignment operator (deleted).
       */
      PolytopeQuadratureIndex& operator=(const PolytopeQuadratureIndex&) = delete;

      /**
       * @brief Move constructor.
       */
      PolytopeQuadratureIndex(PolytopeQuadratureIndex&& other) noexcept
        : m_dimensions(std::move(other.m_dimensions))
      {}

      /**
       * @brief Move assignment operator.
       */
      PolytopeQuadratureIndex& operator=(PolytopeQuadratureIndex&& other) noexcept
      {
        m_dimensions = std::move(other.m_dimensions);
        return *this;
      }

      /**
       * @brief Initializes the index for a mesh of given topological dimension.
       * @param[in] meshDim Topological mesh dimension
       *
       * Allocates dimension buckets for dimensions
       * @f$ 0, \dots, \texttt{meshDim} @f$.
       */
      void initialize(size_t meshDim)
      {
        m_dimensions.resize(meshDim + 1);
      }

      /**
       * @brief Resizes storage for polytopes of dimension @p d.
       * @param[in] d Polytope dimension
       * @param[in] count Number of polytope slots to reserve
       *
       * Ensures that the cache can store entries for at least @p count
       * polytopes of dimension @p d.
       *
       * @note The dimension bucket must already have been initialized.
       */
      void resize(size_t d, size_t count)
      {
        assert(d < m_dimensions.size());
        auto& dim = m_dimensions[d];
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count);
      }

      /**
       * @brief Clears all cached quadratures.
       *
       * Releases all cached quadratures in all dimensions. The index remains
       * initialized, but contains no entries afterwards.
       */
      void clear()
      {
        for (auto& dim : m_dimensions)
        {
          std::unique_lock<std::shared_mutex> wr(dim.mutex);
          dim.slots.clear();
        }
      }

      /**
       * @brief Gets or creates a cached quadrature for a given polytope and formula.
       *
       * @tparam Factory Callable returning `std::unique_ptr<PolytopeQuadrature>`
       * @param[in] p Pair `(dimension, index)` identifying the polytope
       * @param[in] count Expected number of polytopes in that dimension
       * @param[in] qf Quadrature formula used as cache key
       * @param[in] factory Nullary factory called only on cache miss
       * @returns Reference to the cached polytope quadrature
       *
       * Lookup proceeds in three stages:
       * - check the hot entry,
       * - scan the bounded owned cache,
       * - lazily construct and insert on miss.
       *
       * On hits in the bounded cache, the entry is promoted to the hot slot.
       * On misses, a new owned entry is inserted into the bounded cache, using
       * round-robin replacement once the slot is full, and is then promoted to
       * the hot slot.
       *
       * @throws Rodin::Alert::Exception if:
       * - the dimension is outside the initialized range,
       * - the polytope index exceeds the supplied count.
       */
      template <class Factory>
      const PolytopeQuadrature& get(
          const std::pair<size_t, Index>& p,
          size_t count,
          const QF::QuadratureFormulaBase& qf,
          const Factory& factory) const
      {
        const size_t d = p.first;
        const Index idx = p.second;

        if (m_dimensions.size() <= d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Requested dimension " << d
            << " is out of range for an index initialized with "
            << m_dimensions.size() << " dimension buckets."
            << Alert::Raise;
        }

        if (idx >= count)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Polytope index " << idx
            << " is out of range for dimension " << d
            << " with expected count " << count << "."
            << Alert::Raise;
        }

        auto& dim = m_dimensions[d];

        // --------------------------------------------------------------------
        // Fast path: shared lookup in an existing slot
        // --------------------------------------------------------------------
        {
          std::shared_lock<std::shared_mutex> rd(dim.mutex);
          if (idx < dim.slots.size())
          {
            const auto& slot = dim.slots[idx];
            std::shared_lock<std::shared_mutex> rdslot(slot.mutex);

            // Hot entry: fastest repeated-hit path.
            if (slot.hot.qf == &qf)
            {
              assert(slot.hot.ptr);
              return *slot.hot.ptr;
            }

            // Bounded cache lookup.
            for (size_t i = 0; i < slot.size; ++i)
            {
              const auto& entry = slot.entries[i];
              if (entry.valid && entry.qf == &qf)
              {
                assert(entry.owner);
                return *entry.owner;
              }
            }
          }
        }

        // --------------------------------------------------------------------
        // Ensure slot exists
        // --------------------------------------------------------------------
        {
          std::unique_lock<std::shared_mutex> wr(dim.mutex);
          if (dim.slots.size() < count)
            dim.slots.resize(count);
        }

        assert(idx < dim.slots.size());
        auto& slot = dim.slots[idx];

        // --------------------------------------------------------------------
        // Slow path: exclusive slot access
        // --------------------------------------------------------------------
        {
          std::unique_lock<std::shared_mutex> wrslot(slot.mutex);

          // Re-check hot entry under exclusive access.
          if (slot.hot.qf == &qf)
          {
            assert(slot.hot.ptr);
            return *slot.hot.ptr;
          }

          // Re-check bounded cache under exclusive access.
          for (size_t i = 0; i < slot.size; ++i)
          {
            auto& entry = slot.entries[i];
            if (entry.valid && entry.qf == &qf)
            {
              assert(entry.owner);
              slot.promote(entry);
              return *entry.owner;
            }
          }

          // Miss: build and insert into bounded cache.
          CacheEntry& entry = slot.insert(&qf, factory());
          assert(entry.owner);

          // Promote inserted entry to the hot path.
          slot.promote(entry);
          return *entry.owner;
        }
      }

      /**
       * @brief Gets the number of initialized dimension buckets.
       * @returns Number of dimension buckets
       */
      size_t dimensions() const
      {
        return m_dimensions.size();
      }

    private:
      /**
       * @brief Single owned cache entry.
       *
       * This is the owning storage unit for a cached polytope quadrature.
       */
      struct CacheEntry
      {
        const QF::QuadratureFormulaBase* qf = nullptr; ///< Cache key
        std::unique_ptr<PolytopeQuadrature> owner;     ///< Owned cached quadrature
        bool valid = false;                            ///< Whether this entry is populated
      };

      /**
       * @brief Non-owning hot entry.
       *
       * The hot entry mirrors one of the owned cache entries to accelerate the
       * most common repeated-hit case.
       */
      struct HotEntry
      {
        const QF::QuadratureFormulaBase* qf = nullptr; ///< Hot key
        PolytopeQuadrature* ptr = nullptr;             ///< Non-owning pointer to owned entry
      };

      /**
       * @brief Cache slot associated with one concrete polytope.
       *
       * Each slot stores:
       * - one hot entry,
       * - one bounded ring buffer of owned entries.
       *
       * Ownership remains solely in the bounded cache entries. The hot entry is
       * only a mirror of one owned entry.
       */
      struct Slot
      {
        /**
         * @brief Maximum number of owned quadratures cached per polytope.
         *
         * Lookup is linear in this value, so it is intentionally kept small.
         */
        static constexpr size_t Capacity = 8;

        Slot() = default;
        Slot(const Slot&) = delete;
        Slot& operator=(const Slot&) = delete;

        /**
         * @brief Move constructor.
         *
         * The mutex is intentionally default-constructed in the moved instance.
         */
        Slot(Slot&& other) noexcept
          : hot(other.hot),
            entries(std::move(other.entries)),
            size(std::exchange(other.size, 0)),
            next(std::exchange(other.next, 0))
        {}

        /**
         * @brief Move assignment operator.
         *
         * The mutex is intentionally default-constructed in the moved instance.
         */
        Slot& operator=(Slot&& other) noexcept
        {
          hot = other.hot;
          entries = std::move(other.entries);
          size = std::exchange(other.size, 0);
          next = std::exchange(other.next, 0);
          return *this;
        }

        /**
         * @brief Promotes an owned entry to the hot path.
         * @param[in] entry Cache entry to mirror in the hot slot
         *
         * The hot entry never owns. It only mirrors the provided owned entry.
         */
        void promote(const CacheEntry& entry)
        {
          assert(entry.valid);
          assert(entry.qf);
          assert(entry.owner);
          hot.qf = entry.qf;
          hot.ptr = entry.owner.get();
        }

        /**
         * @brief Inserts an entry into the bounded cache.
         * @param[in] qf Quadrature formula key
         * @param[in] owner Newly constructed cached quadrature
         * @returns Reference to the inserted owned entry
         *
         * If the cache is not full, the entry is appended.
         * Otherwise, the next entry is overwritten in round-robin order.
         *
         * If the overwritten entry is currently mirrored by the hot slot,
         * the hot slot is first invalidated to avoid a stale non-owning pointer.
         */
        CacheEntry& insert(
            const QF::QuadratureFormulaBase* qf,
            std::unique_ptr<PolytopeQuadrature> owner)
        {
          assert(qf);
          assert(owner);

          size_t i;
          if (size < Capacity)
          {
            i = size++;
          }
          else
          {
            i = next;
            next = (next + 1) % Capacity;

            // Refinement:
            // If the entry being overwritten is currently mirrored in the hot
            // slot, invalidate the hot entry first to avoid keeping a stale
            // non-owning pointer.
            if (hot.ptr == entries[i].owner.get())
            {
              hot.qf = nullptr;
              hot.ptr = nullptr;
            }
          }

          auto& entry = entries[i];
          entry.qf = qf;
          entry.owner = std::move(owner);
          entry.valid = true;
          return entry;
        }

        HotEntry hot;                                ///< Non-owning fast path
        std::array<CacheEntry, Capacity> entries;    ///< Owned bounded cache
        size_t size = 0;                             ///< Number of valid entries
        size_t next = 0;                             ///< Round-robin eviction pointer
        mutable std::shared_mutex mutex;             ///< Slot-level synchronization
      };

      /**
       * @brief Storage bucket for all polytopes in a fixed dimension.
       *
       * A dimension bucket stores all slots for one polytope dimension.
       */
      struct Dimension
      {
        Dimension() = default;
        Dimension(const Dimension&) = delete;
        Dimension& operator=(const Dimension&) = delete;

        /**
         * @brief Move constructor.
         *
         * The mutex is intentionally default-constructed in the moved instance.
         */
        Dimension(Dimension&& other) noexcept
          : slots(std::move(other.slots))
        {}

        /**
         * @brief Move assignment operator.
         *
         * The mutex is intentionally default-constructed in the moved instance.
         */
        Dimension& operator=(Dimension&& other) noexcept
        {
          slots = std::move(other.slots);
          return *this;
        }

        std::deque<Slot> slots;              ///< Polytope slots for this dimension
        mutable std::shared_mutex mutex;     ///< Dimension-level synchronization
      };

      mutable std::vector<Dimension> m_dimensions; ///< Dimension-partitioned storage
  };
}

#endif
