/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPETRANSFORMATIONINDEX_H
#define RODIN_GEOMETRY_POLYTOPETRANSFORMATIONINDEX_H

/**
 * @file
 * @brief Index for storing polytope transformations.
 */

#include <atomic>
#include <deque>
#include <mutex>
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <shared_mutex>

#include <boost/serialization/access.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include "Rodin/Types.h"
#include "Rodin/Geometry/PolytopeTransformation.h"

namespace Rodin::Geometry
{
  /**
   * @brief Thread-safe index for managing polytope transformations.
   *
   * This class maintains a mapping from polytopes (identified by dimension
   * and index) to their associated geometric transformations. It supports
   * lazy initialization of transformations via factory functions and
   * concurrent access through atomic operations and shared mutexes.
   *
   * The class is designed for high-performance access patterns where
   * transformations are frequently read but infrequently written.
   *
   * @note This class supports Boost serialization for saving and loading
   * transformation data.
   */
  class PolytopeTransformationIndex
  {
    friend class boost::serialization::access;

    /**
     * @brief Storage slot for a single transformation.
     *
     * Contains both an atomic pointer for fast lock-free reads and a
     * unique_ptr for ownership. The atomic pointer is kept synchronized
     * with the unique_ptr.
     */
    struct Slot
    {
      std::atomic<PolytopeTransformation*> ptr{nullptr}; ///< Atomic pointer for fast access
      std::unique_ptr<PolytopeTransformation> owner;     ///< Unique pointer owning the transformation

      /**
       * @brief Serialization save method.
       * @param[in,out] ar Archive object
       * @param[in] version Serialization version (unused)
       */
      template <class Archive>
      void save(Archive& ar, const unsigned int) const
      {
        ar & owner; // polymorphic unique_ptr
      }

      /**
       * @brief Serialization load method.
       * @param[in,out] ar Archive object
       * @param[in] version Serialization version (unused)
       *
       * Restores both the unique_ptr and synchronizes the atomic pointer.
       */
      template <class Archive>
      void load(Archive& ar, const unsigned int)
      {
        ar & owner;
        ptr.store(owner.get(), std::memory_order_relaxed);
      }

      BOOST_SERIALIZATION_SPLIT_MEMBER()
    };

    /**
     * @brief Storage for transformations of polytopes in a single dimension.
     *
     * Contains a deque of transformation slots along with a shared mutex
     * for thread-safe concurrent access.
     */
    struct Dimension
    {
      std::deque<Slot> slots;          ///< Storage for transformation slots
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
        : slots(std::move(other.slots)) {}

      /**
       * @brief Move assignment operator.
       */
      Dimension& operator=(Dimension&& other) noexcept
      {
        slots = std::move(other.slots);
        return *this;
      }

      /**
       * @brief Serialization method.
       * @param[in,out] ar Archive object
       * @param[in] version Serialization version (unused)
       *
       * @note The mutex is not serialized.
       */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & slots; // mutex is not serialized
      }
    };

  public:
    /**
     * @brief Default constructor.
     */
    PolytopeTransformationIndex() = default;

    /**
     * @brief Destructor.
     */
    ~PolytopeTransformationIndex() = default;

    /**
     * @brief Copy constructor (deleted).
     */
    PolytopeTransformationIndex(const PolytopeTransformationIndex&) = delete;

    /**
     * @brief Copy assignment operator (deleted).
     */
    PolytopeTransformationIndex& operator=(const PolytopeTransformationIndex&) = delete;

    /**
     * @brief Move constructor.
     */
    PolytopeTransformationIndex(PolytopeTransformationIndex&& other) noexcept
      : m_dimensions(std::move(other.m_dimensions))
    {}

    /**
     * @brief Move assignment operator.
     */
    PolytopeTransformationIndex& operator=(PolytopeTransformationIndex&& other) noexcept
    {
      m_dimensions = std::move(other.m_dimensions);
      return *this;
    }

    /**
     * @brief Initializes the index for a mesh of given dimension.
     * @param[in] meshDim Topological dimension of the mesh
     *
     * Must be called before any other operations. Allocates storage for
     * dimensions 0 through @p meshDim.
     */
    void initialize(size_t meshDim)
    {
      m_dimensions.resize(meshDim + 1);
    }

    /**
     * @brief Resizes storage for polytopes of dimension @p d.
     * @param[in] d Dimension of polytopes
     * @param[in] count Number of polytopes to allocate space for
     *
     * Ensures the internal storage can hold at least @p count polytopes
     * of dimension @p d.
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
     * @brief Sets the transformation for a polytope.
     * @param[in] p Pair of (dimension, index) identifying the polytope
     * @param[in] obj Unique pointer to the transformation
     *
     * Stores the transformation and updates the atomic pointer for fast access.
     * Automatically resizes storage if needed.
     */
    void set(const std::pair<size_t, Index>& p,
             std::unique_ptr<PolytopeTransformation> obj)
    {
      const size_t d   = p.first;
      const Index  idx = p.second;

      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];
      std::unique_lock<std::shared_mutex> wr(dim.mutex);
      if (dim.slots.size() <= idx)
        dim.slots.resize(idx + 1);

      assert(idx < dim.slots.size());
      Slot& s = dim.slots[idx];
      s.owner = std::move(obj);
      s.ptr.store(s.owner.get(), std::memory_order_release);
    }

    /**
     * @brief Sets the transformation for a polytope with known count.
     * @param[in] p Pair of (dimension, index) identifying the polytope
     * @param[in] count Expected number of polytopes in this dimension
     * @param[in] obj Unique pointer to the transformation
     *
     * More efficient than set() when the total count is known in advance.
     */
    void set(const std::pair<size_t, Index>& p,
             size_t count,
             std::unique_ptr<PolytopeTransformation> obj)
    {
      const size_t d   = p.first;
      const Index  idx = p.second;

      assert(idx < count);

      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];
      std::unique_lock<std::shared_mutex> wr(dim.mutex);
      if (dim.slots.size() < count)
        dim.slots.resize(count);

      assert(idx < dim.slots.size());
      Slot& s = dim.slots[idx];
      s.owner = std::move(obj);
      s.ptr.store(s.owner.get(), std::memory_order_release);
    }

    /**
     * @brief Gets or creates a transformation using a factory.
     * @tparam Factory Callable type that creates transformations
     * @param[in] p Pair of (dimension, index) identifying the polytope
     * @param[in] count Expected number of polytopes in this dimension
     * @param[in] factory Factory function to create transformation if needed
     * @returns Reference to the transformation
     *
     * Uses a two-phase locking strategy: first tries a fast shared-lock read,
     * then upgrades to exclusive lock only if the transformation needs to be
     * created. The factory is called at most once per polytope.
     *
     * @note The Factory must be callable with signature:
     *       `std::unique_ptr<PolytopeTransformation>(size_t d, Index idx)`
     */
    template <class Factory>
    const PolytopeTransformation&
    get(const std::pair<size_t, Index>& p, size_t count, const Factory& factory) const
    {
      const size_t d   = p.first;
      const Index  idx = p.second;

      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];

      // Fast path: read under shared lock
      {
        std::shared_lock<std::shared_mutex> rd(dim.mutex);
        if (idx < dim.slots.size())
        {
          assert(idx < count);
          const Slot& s = dim.slots[idx];
          if (auto* q = s.ptr.load(std::memory_order_acquire)) return *q;
        }
      }

      // Slow path: exclusive lock, ensure size, init in-place
      {
        std::unique_lock<std::shared_mutex> wr(dim.mutex);

        if (dim.slots.size() < count)
          dim.slots.resize(count);

        assert(idx < dim.slots.size());
        Slot& s = dim.slots[idx];
        PolytopeTransformation* q = s.ptr.load(std::memory_order_relaxed);
        if (!q)
        {
          auto up = factory(d, idx);
          s.owner = std::move(up);
          q = s.owner.get();
          s.ptr.store(q, std::memory_order_release);
        }
        return *q; // still under wr; safe
      }
    }

    /**
     * @brief Clears all stored transformations.
     *
     * Releases all transformation objects and resets internal storage.
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

    /**
     * @brief Serialization save method.
     * @param[in,out] ar Archive object
     * @param[in] version Serialization version (unused)
     */
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
      ar & m_dimensions;
    }

    /**
     * @brief Serialization load method.
     * @param[in,out] ar Archive object
     * @param[in] version Serialization version (unused)
     *
     * Clears existing data before loading.
     */
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
      clear();
      ar & m_dimensions;
      // ptr fields restored by Slot::load
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

  private:
    mutable std::vector<Dimension> m_dimensions; ///< Storage for each dimension
  };
} // namespace Rodin::Geometry

#endif // RODIN_GEOMETRY_POLYTOPETRANSFORMATIONINDEX_H
