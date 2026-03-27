/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_QUADRATUREPOINTINDEX_H
#define RODIN_GEOMETRY_QUADRATUREPOINTINDEX_H

/**
 * @file
 * @brief Thread-safe lazy index for caching quadrature points per polytope.
 */

#include <atomic>
#include <deque>
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <shared_mutex>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Geometry
{
  /**
   * @brief Thread-safe index for caching quadrature points per polytope.
   *
   * This class maintains a mapping from polytopes (identified by dimension
   * and index) to their associated quadrature points for a given quadrature
   * order and geometry type. It supports lazy initialization of point sets
   * via factory functions and concurrent access through atomic operations
   * and shared mutexes.
   *
   * Each cached entry stores a vector of Geometry::Point objects with
   * pre-populated geometric caches (Jacobian, distortion, etc.) so that
   * subsequent const reads are safe for concurrent access without mutation.
   *
   * @note Not copyable. Move-only.
   * @see PolytopeTransformationIndex, AttributeIndex
   */
  class QuadraturePointIndex
  {
    /**
     * @brief A set of quadrature points cached for one polytope.
     *
     * Stores the quadrature order and geometry type used to create the
     * points, so that stale entries can be detected when the requested
     * order or geometry changes.
     */
    struct PointSet
    {
      std::vector<Point> points;  ///< Cached quadrature points
      size_t order;               ///< Quadrature order used
      Polytope::Type geometry;    ///< Polytope geometry type
    };

    /**
     * @brief Storage slot for a single polytope's quadrature points.
     *
     * Contains both an atomic pointer for fast lock-free reads and a
     * unique_ptr for ownership. The atomic pointer is kept synchronized
     * with the unique_ptr.
     */
    struct Slot
    {
      std::atomic<PointSet*> ptr{nullptr}; ///< Atomic pointer for fast access
      std::unique_ptr<PointSet> owner;     ///< Unique pointer owning the PointSet
    };

    /**
     * @brief Storage for quadrature points of polytopes in a single dimension.
     *
     * Contains a deque of slots along with a shared mutex for thread-safe
     * concurrent access.
     */
    struct Dimension
    {
      std::deque<Slot> slots;          ///< Storage for quadrature point slots
      mutable std::shared_mutex mutex; ///< Mutex for thread-safe access

      /**
       * @brief Default constructor.
       */
      Dimension() = default;

      /**
       * @brief Copy constructor (deleted).
       */
      Dimension(const Dimension&) = delete;

      /**
       * @brief Copy assignment operator (deleted).
       */
      Dimension& operator=(const Dimension&) = delete;

      /**
       * @brief Move constructor.
       */
      Dimension(Dimension&& other) noexcept
        : slots(std::move(other.slots))
      {}

      /**
       * @brief Move assignment operator.
       */
      Dimension& operator=(Dimension&& other) noexcept
      {
        slots = std::move(other.slots);
        return *this;
      }
    };

  public:
    /**
     * @brief Default constructor.
     */
    QuadraturePointIndex() = default;

    /**
     * @brief Destructor.
     */
    ~QuadraturePointIndex() = default;

    /**
     * @brief Copy constructor (deleted).
     */
    QuadraturePointIndex(const QuadraturePointIndex&) = delete;

    /**
     * @brief Copy assignment operator (deleted).
     */
    QuadraturePointIndex& operator=(const QuadraturePointIndex&) = delete;

    /**
     * @brief Move constructor.
     */
    QuadraturePointIndex(QuadraturePointIndex&& other) noexcept
      : m_dimensions(std::move(other.m_dimensions))
    {}

    /**
     * @brief Move assignment operator.
     */
    QuadraturePointIndex& operator=(QuadraturePointIndex&& other) noexcept
    {
      m_dimensions = std::move(other.m_dimensions);
      return *this;
    }

    /**
     * @brief Initializes the index for a mesh of given dimension.
     * @param[in] meshDim Topological dimension of the mesh
     *
     * Allocates storage for dimensions 0 through @p meshDim.
     * Idempotent: if already initialized with at least this dimension, no-op.
     */
    void initialize(size_t meshDim)
    {
      const size_t newSize = meshDim + 1;
      if (m_dimensions.size() >= newSize)
        return;
      m_dimensions.resize(newSize);
    }

    /**
     * @brief Resizes storage for polytopes of dimension @p d.
     * @param[in] d Dimension of polytopes
     * @param[in] count Number of polytopes to allocate space for
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
     * @brief Gets or creates cached quadrature points for a polytope.
     * @tparam Factory Callable type that creates quadrature points
     * @param[in] d Polytope dimension
     * @param[in] idx Polytope index
     * @param[in] order Quadrature order
     * @param[in] geometry Polytope geometry type
     * @param[in] count Expected number of polytopes in this dimension
     * @param[in] factory Factory function to create points if needed
     * @returns Const reference to the cached vector of quadrature points
     *
     * Uses a two-phase locking strategy: first tries a fast shared-lock read,
     * then upgrades to exclusive lock only if the points need to be created.
     *
     * If a cached entry exists with matching order and geometry, it is returned
     * immediately. Otherwise, the factory is called to create a new set of
     * points. The factory must pre-populate all Point caches (getDistortion(),
     * etc.) so that subsequent const access is thread-safe.
     *
     * @note The Factory must be callable with signature:
     *       `std::vector<Point>(size_t d, Index idx)`
     */
    template <class Factory>
    const std::vector<Point>&
    get(size_t d, Index idx,
        size_t order, Polytope::Type geometry,
        size_t count, const Factory& factory) const
    {
      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];

      // Fast path: read under shared lock
      {
        std::shared_lock<std::shared_mutex> rd(dim.mutex);
        if (idx < dim.slots.size())
        {
          const Slot& s = dim.slots[idx];
          if (auto* ps = s.ptr.load(std::memory_order_acquire))
          {
            if (ps->order == order && ps->geometry == geometry)
              return ps->points;
          }
        }
      }

      // Slow path: exclusive lock, create/replace
      {
        std::unique_lock<std::shared_mutex> wr(dim.mutex);

        if (dim.slots.size() < count)
          dim.slots.resize(count);

        assert(idx < dim.slots.size());
        Slot& s = dim.slots[idx];
        PointSet* ps = s.ptr.load(std::memory_order_relaxed);
        if (!ps || ps->order != order || ps->geometry != geometry)
        {
          auto up = std::make_unique<PointSet>();
          up->order = order;
          up->geometry = geometry;
          up->points = factory(d, idx);
          s.owner = std::move(up);
          ps = s.owner.get();
          s.ptr.store(ps, std::memory_order_release);
        }
        return ps->points;
      }
    }

    /**
     * @brief Clears all cached quadrature points.
     *
     * Releases all cached data and resets internal storage.
     * Thread-safe with respect to concurrent operations.
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
     * @brief Gets the number of dimensions.
     * @returns Number of polytope dimensions managed by this index
     */
    size_t dimensions() const
    {
      return m_dimensions.size();
    }

  private:
    mutable std::vector<Dimension> m_dimensions; ///< Storage for each dimension
  };
} // namespace Rodin::Geometry

#endif // RODIN_GEOMETRY_QUADRATUREPOINTINDEX_H
