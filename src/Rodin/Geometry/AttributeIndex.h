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
#include <vector>
#include <shared_mutex>

#include <boost/serialization/access.hpp>

#include "Rodin/Geometry/Types.h" // Attribute, Index, FlatSet, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE

namespace Rodin::Geometry
{
  /**
   * @brief Manages attribute assignments for mesh polytopes.
   *
   * This class maintains a mapping from polytopes (identified by dimension
   * and index) to their associated attributes. Attributes are typically used
   * to identify material regions, boundary markers, or other domain-specific
   * properties.
   *
   * The class is thread-safe for concurrent access through the use of
   * shared mutexes for each dimension.
   *
   * @note This class supports Boost serialization for saving and loading
   * mesh attribute information.
   */
  class AttributeIndex
  {
    friend class boost::serialization::access;

    /**
     * @brief Storage for attributes of polytopes in a single dimension.
     *
     * Contains a vector of attributes indexed by polytope index, along with
     * a shared mutex for thread-safe access.
     */
    struct Dimension
    {
      friend class boost::serialization::access;

      std::vector<Attribute> slots; ///< Attribute values indexed by polytope index

      mutable std::shared_mutex mutex; ///< Mutex for thread-safe access

      /**
       * @brief Default constructor.
       */
      Dimension() = default;

      /**
       * @brief Copy constructor.
       */
      Dimension(const Dimension& other)
        : slots(other.slots)
      {}

      /**
       * @brief Copy assignment operator.
       */
      Dimension& operator=(const Dimension& other)
      {
        if (this != &other)
        {
          slots = other.slots;
        }
        return *this;
      }

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
       */
      AttributeIndex() = default;

      /**
       * @brief Destructor.
       */
      ~AttributeIndex() = default;

      /**
       * @brief Copy constructor.
       */
      AttributeIndex(const AttributeIndex& other)
        : m_dimensions(other.m_dimensions)
      {}

      /**
       * @brief Copy assignment operator (deleted).
       */
      AttributeIndex& operator=(const AttributeIndex&) = delete;

      /**
       * @brief Move constructor.
       */
      AttributeIndex(AttributeIndex&& other) noexcept
        : m_dimensions(std::move(other.m_dimensions))
      {}

      /**
       * @brief Move assignment operator.
       */
      AttributeIndex& operator=(AttributeIndex&& other) noexcept
      {
        m_dimensions = std::move(other.m_dimensions);
        return *this;
      }

      /**
       * @brief Initializes the attribute index for a mesh of given dimension.
       * @param[in] meshDim Topological dimension of the mesh
       *
       * Must be called once before any concurrent access. Allocates storage
       * for dimensions 0 through @p meshDim.
       */
      void initialize(size_t meshDim)
      {
        m_dimensions.resize(meshDim + 1);
      }

      /**
       * @brief Resizes the storage for polytopes of dimension @p d.
       * @param[in] d Dimension of polytopes
       * @param[in] count Number of polytopes to allocate space for
       *
       * Ensures that the internal storage can hold at least @p count
       * polytopes of dimension @p d.
       */
      void resize(size_t d, size_t count)
      {
        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
      }

      /**
       * @brief Sets the attribute for a polytope.
       * @param[in] p Pair of (dimension, index) identifying the polytope
       * @param[in] attr Attribute value to assign
       *
       * Automatically resizes storage if needed to accommodate the polytope.
       */
      void set(const std::pair<size_t, Index>& p, Attribute attr)
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());

        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() <= idx)
          dim.slots.resize(idx + 1, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
        assert(idx < dim.slots.size());
        dim.slots[idx] = attr;
      }

      /**
       * @brief Sets the attribute for a polytope with known count.
       * @param[in] p Pair of (dimension, index) identifying the polytope
       * @param[in] count Expected number of polytopes in this dimension
       * @param[in] attr Attribute value to assign
       *
       * More efficient than set() when the total count is known in advance.
       */
      void set(const std::pair<size_t, Index>& p, size_t count, Attribute attr)
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());
        assert(idx < count);

        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
        dim.slots[idx] = attr;
      }

      /**
       * @brief Gets the attribute for a polytope.
       * @param[in] p Pair of (dimension, index) identifying the polytope
       * @param[in] count Expected number of polytopes in this dimension
       * @returns The attribute value
       *
       * Returns the current value, growing storage to @p count if needed.
       * Uses shared locking for read access, upgrading to exclusive locking
       * only if resizing is required.
       */
      Attribute get(const std::pair<size_t, Index>& p, size_t count) const
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());
        assert(idx < count);

        // Read phase
        {
          const Dimension& dim = m_dimensions[d];                  // const ref
          std::shared_lock<std::shared_mutex> rd(dim.mutex);
          if (idx < dim.slots.size())
            return dim.slots[idx];
        }

        // Write phase
        {
          Dimension& dim = m_dimensions[d];                        // non-const ref
          std::unique_lock<std::shared_mutex> wr(dim.mutex);
          if (dim.slots.size() < count)
            dim.slots.resize(count, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
          if (idx >= dim.slots.size())
            dim.slots.resize(idx + 1, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
          return dim.slots[idx];
        }
      }

      /**
       * @brief Gets all unique attributes in a given dimension.
       * @param[in] d Dimension of polytopes
       * @returns Set of all unique attribute values
       *
       * Useful for identifying all material regions or markers present
       * in the mesh.
       */
      FlatSet<Attribute> getAttributes(size_t d) const
      {
        FlatSet<Attribute> out;
        const auto& dim = m_dimensions.at(d);
        std::shared_lock<std::shared_mutex> rd(dim.mutex);
        for (Attribute a : dim.slots)
          out.insert(a);
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
      mutable std::vector<Dimension> m_dimensions; ///< Storage for each dimension
  };
}
#endif
