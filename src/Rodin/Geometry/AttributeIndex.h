/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ATTRIBUTEINDEX_H
#define RODIN_GEOMETRY_ATTRIBUTEINDEX_H

/**
 * @file
 * @brief Attribute indexing for mesh polytopes.
 */

#include <mutex>
#include <optional>
#include <vector>
#include <shared_mutex>

#include <boost/serialization/access.hpp>

#include "Rodin/Geometry/Types.h" // Attribute, Index, FlatSet

namespace Rodin::Geometry
{
  /**
   * @brief Attribute indexing for mesh polytopes.
   *
   * This class stores, for each topological dimension @c d, an array of optional
   * attributes indexed by polytope index. Attributes are typically used as
   * material-region ids, boundary markers, or any other user-defined tags.
   *
   * Data layout
   * -----------
   * The storage is organized by dimension:
   * - For each @c d in @c [0, meshDim], a @ref Dimension stores:
   *   - @c slots: a vector of @c Optional<Attribute>, indexed by polytope index.
   *   - @c mutex: a @c std::shared_mutex protecting accesses to @c slots.
   *
   * Thread-Safety Model
   * -------------------
   * This class provides fine-grained thread-safety at the level of individual
   * dimensions:
   * - Each dimension owns an independent @c std::shared_mutex.
   * - Operations that read/write attributes in a given dimension synchronize
   *   using that dimension's mutex.
   *
   * After initialization, concurrent calls to the following methods are safe
   * w.r.t. each other (subject to each method's preconditions):
   * - @ref resize(size_t,size_t)
   * - @ref set(const std::pair<size_t,Index>&,Attribute)
   * - @ref set(const std::pair<size_t,Index>&,size_t,Attribute)
   * - @ref unset(const std::pair<size_t,Index>&,size_t)
   * - @ref get(const std::pair<size_t,Index>&,size_t) const
   * - @ref getAttributes(size_t) const
   *
   * Structural Immutability
   * -----------------------
   * After @ref initialize() completes, the number of stored dimensions
   * (i.e. @c m_dimensions.size()) is immutable for the lifetime of the object.
   * This avoids requiring a global lock on the dimension vector itself.
   *
   * Therefore:
   * - @ref initialize() must be called before any concurrent access.
   * - @ref initialize() is not thread-safe and must not be called concurrently.
   * - The mesh dimension must remain constant for the lifetime of the object.
   *
   * Copy and Move Semantics
   * -----------------------
   * Copy construction:
   * - Copies per-dimension storage using the source per-dimension locks.
   *
   * Copy assignment:
   * - Requires both objects to be initialized with the same number of dimensions.
   * - Copies per dimension using the per-dimension locks.
   * - Is NOT an atomic snapshot across all dimensions: a concurrent reader may
   *   observe a mix of old and new dimension contents while assignment is in
   *   progress.
   *
   * Move construction / move assignment:
   * - Transfer ownership of the internal dimension vector.
   * - Require external synchronization: no other thread may access either object
   *   concurrently with the move.
   *
   * Usage Contract
   * --------------
   * 1. Call @ref initialize(meshDim) exactly once (or idempotently with the same
   *    value) before any concurrent use.
   * 2. Do not change the mesh dimension after initialization.
   * 3. Ensure storage exists (via @ref resize() or via @ref set()) before relying
   *    on reads returning non-null values.
   * 4. Provide external synchronization around move operations.
   */
  class AttributeIndex
  {
    friend class boost::serialization::access;

    /**
     * @brief Storage for attributes of polytopes in a single dimension.
     *
     * Contains:
     * - @c slots: attribute values indexed by polytope index.
     * - @c mutex: a shared mutex protecting @c slots.
     *
     * Thread-safety:
     * - Multiple concurrent readers are allowed via shared locking.
     * - Writers use exclusive locking.
     */
    struct Dimension
    {
      friend class boost::serialization::access;

      std::vector<Optional<Attribute>> slots; ///< Attribute values indexed by polytope index

      mutable std::shared_mutex mutex; ///< Mutex for thread-safe access

      /**
       * @brief Default constructor.
       */
      Dimension() = default;

      /**
       * @brief Copy constructor.
       *
       * Copies the slot vector from @p other under @p other's shared lock.
       *
       * Thread-safety:
       * - Safe w.r.t. concurrent readers/writers of @p other.
       */
      Dimension(const Dimension& other)
      {
        std::shared_lock<std::shared_mutex> rd(other.mutex);
        slots = other.slots;
      }

      /**
       * @brief Copy assignment operator.
       *
       * Assigns slots from @p other to @c *this using per-object locking.
       * A deadlock-avoiding lock ordering is used based on object addresses.
       *
       * Thread-safety:
       * - Safe w.r.t. concurrent readers/writers of either dimension.
       * - Not an atomic snapshot w.r.t. external observers of the whole container.
       */
      Dimension& operator=(const Dimension& other)
      {
        if (this == &other)
          return *this;

        // Lock ordering by address to avoid deadlock.
        const Dimension* first  = this;
        const Dimension* second = &other;
        if (std::less<const Dimension*>{}(second, first))
          std::swap(first, second);

        if (first == this)
        {
          std::unique_lock<std::shared_mutex> lock_this(mutex);
          std::shared_lock<std::shared_mutex> lock_other(other.mutex);
          slots = other.slots;
        }
        else
        {
          std::shared_lock<std::shared_mutex> lock_other(other.mutex);
          std::unique_lock<std::shared_mutex> lock_this(mutex);
          slots = other.slots;
        }

        return *this;
      }

      /**
       * @brief Move constructor.
       *
       * Moves the slot vector from @p other under @p other's exclusive lock.
       *
       * Thread-safety:
       * - Requires that no other thread accesses @p other concurrently.
       */
      Dimension(Dimension&& other) noexcept
      {
        std::unique_lock<std::shared_mutex> wr(other.mutex);
        slots = std::move(other.slots);
      }

      /**
       * @brief Move assignment operator.
       *
       * Moves slots from @p other into @c *this using exclusive locks on both.
       *
       * Thread-safety:
       * - Requires external synchronization: no concurrent access to either
       *   dimension while moving.
       */
      Dimension& operator=(Dimension&& other) noexcept
      {
        if (this == &other)
          return *this;

        std::unique_lock<std::shared_mutex> this_lock(mutex, std::defer_lock);
        std::unique_lock<std::shared_mutex> other_lock(other.mutex, std::defer_lock);
        std::lock(this_lock, other_lock);

        slots = std::move(other.slots);
        return *this;
      }

      /**
       * @brief Serialization method for Boost.Serialization.
       * @param[in,out] ar Archive object
       * @param[in] version Serialization version (unused)
       */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & slots;
      }
    };

    public:
      /**
       * @brief Default constructor.
       *
       * Constructs an uninitialized index. Call @ref initialize() before use.
       */
      AttributeIndex() = default;

      /**
       * @brief Destructor.
       */
      ~AttributeIndex() = default;

      /**
       * @brief Copy constructor.
       *
       * Copies all dimensions from @p other.
       *
       * Thread-safety:
       * - Safe if @p other is not being moved concurrently.
       * - Copy is not a globally atomic snapshot across dimensions.
       */
      AttributeIndex(const AttributeIndex& other)
        : m_dimensions(other.m_dimensions)
      {}

      /**
       * @brief Copy assignment operator.
       *
       * Copies all dimensions from @p other into @c *this.
       *
       * Preconditions:
       * - Both objects must already be initialized with the same mesh dimension
       *   (i.e. same @c m_dimensions.size()).
       *
       * Thread-safety:
       * - Uses per-dimension locks.
       * - Not an atomic snapshot across dimensions.
       */
      AttributeIndex& operator=(const AttributeIndex& other)
      {
        if (this == &other)
          return *this;

        // Contract: both already initialized with the same mesh dimension.
        assert(m_dimensions.size() == other.m_dimensions.size());

        for (size_t d = 0; d < m_dimensions.size(); ++d)
          m_dimensions[d] = other.m_dimensions[d]; // Dimension::operator= locks

        return *this;
      }

      /**
       * @brief Move constructor.
       *
       * Transfers internal dimension storage from @p other.
       *
       * Thread-safety:
       * - Requires external synchronization: no concurrent access to either object.
       */
      AttributeIndex(AttributeIndex&& other) noexcept
        : m_dimensions(std::move(other.m_dimensions))
      {}

      /**
       * @brief Move assignment operator.
       *
       * Transfers internal dimension storage from @p other into @c *this.
       *
       * Thread-safety:
       * - Requires external synchronization: no concurrent access to either object.
       */
      AttributeIndex& operator=(AttributeIndex&& other) noexcept
      {
        if (this == &other)
          return *this;
        m_dimensions = std::move(other.m_dimensions);
        return *this;
      }

      /**
       * @brief Initializes the index for a mesh of topological dimension @p meshDim.
       *
       * Allocates per-dimension storage for dimensions @c 0..meshDim.
       *
       * Structural immutability:
       * - After initialization, the dimension count is immutable.
       * - Calling initialize again is only valid if called idempotently with the
       *   same @p meshDim.
       *
       * Thread-safety:
       * - Not thread-safe. Must be called before any concurrent access.
       *
       * @param[in] meshDim Topological dimension of the mesh.
       */
      void initialize(size_t meshDim)
      {
        const size_t newSize = meshDim + 1;

        // Contract enforcement: after initialization, dimension count is immutable.
        if (!m_dimensions.empty())
        {
          assert(m_dimensions.size() == newSize);
          return;
        }

        m_dimensions.resize(newSize);
      }

      /**
       * @brief Ensures storage capacity for dimension @p d.
       *
       * Grows the slot vector for dimension @p d to at least @p count, filling new
       * entries with @c std::nullopt.
       *
       * Thread-safety:
       * - Thread-safe w.r.t. other operations on the same dimension (protected by
       *   that dimension's mutex).
       *
       * @param[in] d Dimension of polytopes.
       * @param[in] count Minimum number of slots to allocate.
       */
      void resize(size_t d, size_t count)
      {
        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, std::nullopt);
      }

      /**
       * @brief Sets the attribute for a polytope.
       *
       * If storage is insufficient, grows the slot vector to include @p idx and
       * assigns @p attr.
       *
       * Thread-safety:
       * - Thread-safe (dimension-level exclusive lock).
       *
       * @param[in] p Pair (dimension, index) identifying the polytope.
       * @param[in] attr Attribute value to assign.
       */
      void set(const std::pair<size_t, Index>& p, const Optional<Attribute>& attr)
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());

        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() <= idx)
          dim.slots.resize(idx + 1, std::nullopt);
        assert(idx < dim.slots.size());
        dim.slots[idx] = attr;
      }

      /**
       * @brief Sets the attribute for a polytope assuming a known count.
       *
       * Ensures the dimension has at least @p count slots (growing if necessary),
       * then assigns @p attr at @p idx.
       *
       * Preconditions:
       * - @p idx < @p count.
       *
       * Thread-safety:
       * - Thread-safe (dimension-level exclusive lock).
       *
       * @param[in] p Pair (dimension, index) identifying the polytope.
       * @param[in] count Expected number of polytopes in this dimension.
       * @param[in] attr Attribute value to assign.
       */
      void set(const std::pair<size_t, Index>& p, size_t count, const Optional<Attribute>& attr)
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());
        assert(idx < count);

        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, std::nullopt);
        dim.slots[idx] = attr;
      }

      /**
       * @brief Clears (unsets) the attribute for a polytope.
       *
       * Ensures the dimension has at least @p count slots (growing if necessary),
       * then resets the entry at @p idx to @c std::nullopt.
       *
       * Thread-safety:
       * - Thread-safe (dimension-level exclusive lock).
       *
       * @param[in] p Pair (dimension, index) identifying the polytope.
       * @param[in] count Expected number of polytopes in this dimension.
       */
      void unset(const std::pair<size_t, Index>& p, size_t count)
      {
        const auto& [d, idx] = p;
        auto& dim = m_dimensions.at(d);
        assert(d < m_dimensions.size());
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, std::nullopt);
        dim.slots[idx].reset();
      }

      /**
       * @brief Gets the attribute for a polytope.
       *
       * Reads the slot at @p idx under a shared lock.
       *
       * Notes:
       * - This method does not resize. If @p idx is out of range, returns
       *   @c std::nullopt.
       *
       * Preconditions:
       * - @p idx < @p count (debug-checked).
       *
       * Thread-safety:
       * - Thread-safe (dimension-level shared lock).
       *
       * @param[in] p Pair (dimension, index) identifying the polytope.
       * @param[in] count Expected number of polytopes in this dimension.
       * @return The optional attribute value.
       */
      Optional<Attribute> get(const std::pair<size_t, Index>& p, size_t count) const
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());
        assert(idx < count);

        const auto& dim = m_dimensions.at(d);
        std::shared_lock<std::shared_mutex> rd(dim.mutex);

        // No resizing here: assume resize(d, count) happened earlier.
        if (idx >= dim.slots.size())
          return std::nullopt;

        return dim.slots[idx];
      }

      /**
       * @brief Returns the set of all attributes present in dimension @p d.
       *
       * Iterates over all slots in dimension @p d and collects all engaged
       * attributes.
       *
       * Thread-safety:
       * - Thread-safe (dimension-level shared lock).
       *
       * @param[in] d Dimension of polytopes.
       * @return A flat set of unique attribute values.
       */
      FlatSet<Attribute> getAttributes(size_t d) const
      {
        FlatSet<Attribute> out;
        const auto& dim = m_dimensions.at(d);
        std::shared_lock<std::shared_mutex> rd(dim.mutex);
        for (const auto& oa : dim.slots)
        {
          if (oa)
            out.insert(*oa);
        }
        return out;
      }

      /**
       * @brief Serialization method for Boost.Serialization.
       * @param[in,out] ar Archive object
       */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & m_dimensions;
      }

    private:
      /// Storage for each dimension (size fixed after initialize()).
      mutable std::vector<Dimension> m_dimensions;
  };
}
#endif
