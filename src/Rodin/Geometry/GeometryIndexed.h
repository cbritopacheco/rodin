/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GEOMETRYINDEXED_H
#define RODIN_VARIATIONAL_GEOMETRYINDEXED_H

/**
 * @file
 * @brief Container for mapping polytope types to associated data.
 */

#include <array>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <boost/serialization/access.hpp>

#include "Polytope.h"

namespace Rodin::Geometry
{
  /**
   * @brief Container that maps polytope types to associated data.
   *
   * This class provides an efficient mapping from polytope geometry types 
   * (Point, Segment, Triangle, etc.) to arbitrary data of type T. It uses
   * compile-time knowledge of the number of polytope types to provide 
   * O(1) access with minimal memory overhead.
   *
   * # Usage Example
   * @code{.cpp}
   * // Store quadrature orders for each polytope type
   * GeometryIndexed<int> quadratureOrders = {
   *   {Polytope::Type::Point, 1},
   *   {Polytope::Type::Segment, 2},
   *   {Polytope::Type::Triangle, 3},
   *   {Polytope::Type::Quadrilateral, 4},
   *   {Polytope::Type::Tetrahedron, 5},
   *   {Polytope::Type::Hexahedron, 6}
   *   {Polytope::Type::Wedge, 7}
   * };
   * int order = quadratureOrders[Polytope::Type::Triangle]; // Returns 3
   * @endcode
   *
   * @tparam T Type of data to associate with each polytope type
   *
   * @see Polytope::Type
   */
  template <class T>
  class GeometryIndexed
  {
    /// Number of different polytope types
    static constexpr size_t Count = Polytope::Types.size();

    /// Storage type for proper alignment
    using Storage = typename std::aligned_storage<sizeof(T), alignof(T)>::type;

    friend class boost::serialization::access;

  public:
    /**
     * @brief Default constructor that initializes all entries with default values.
     */
    GeometryIndexed()
    {
      for (size_t i = 0; i < Count; ++i)
        new (&m_map[i]) T();
    }

    /**
     * @brief Constructs from an initializer list of polytope type-value pairs.
     * @param l Initializer list of (polytope type, value) pairs
     * @note The list must contain exactly one entry for each polytope type
     */
    GeometryIndexed(std::initializer_list<std::pair<Polytope::Type, T>> l)
    {
      assert(l.size() == Count);
      for (const auto& [type, value] : l)
        new (&m_map[static_cast<size_t>(type)]) T(std::move(value));
    }

    /**
     * @brief Copy constructor.
     * @param other GeometryIndexed object to copy from
     */
    GeometryIndexed(const GeometryIndexed& other)
    {
      for (size_t i = 0; i < Count; ++i)
        new (&m_map[i]) T(*other.ptr(i));
    }

    /**
     * @brief Move constructor.
     * @param other GeometryIndexed object to move from
     */
    GeometryIndexed(GeometryIndexed&& other) noexcept(std::is_nothrow_move_constructible<T>::value)
    {
      for (size_t i = 0; i < Count; ++i)
        new (&m_map[i]) T(std::move(*other.ptr(i)));
    }

    /**
     * @brief Destructor that properly destroys all stored objects.
     */
    ~GeometryIndexed()
    {
      for (size_t i = 0; i < Count; ++i)
        ptr(i)->~T();
    }

    /**
     * @brief Copy assignment operator.
     * @param other GeometryIndexed object to copy from
     * @return Reference to this object
     */
    GeometryIndexed& operator=(const GeometryIndexed& other)
      noexcept(std::is_nothrow_copy_assignable<T>::value)
    {
      if (this != &other)
        for (size_t i = 0; i < Count; ++i)
          *ptr(i) = *other.ptr(i);
      return *this;
    }

    /**
     * @brief Move assignment operator.
     * @param other GeometryIndexed object to move from
     * @return Reference to this object
     */
    GeometryIndexed& operator=(GeometryIndexed&& other)
      noexcept(std::is_nothrow_move_assignable<T>::value)
    {
      if (this != &other)
        for (size_t i = 0; i < Count; ++i)
          *ptr(i) = std::move(*other.ptr(i));
      return *this;
    }

    /**
     * @brief Access operator for mutable access to data by polytope type.
     * @param type Polytope type to access
     * @return Reference to the associated data
     */
    T& operator[](Polytope::Type type)
    {
      return *ptr(static_cast<size_t>(type));
    }

    /**
     * @brief Access operator for const access to data by polytope type.
     * @param type Polytope type to access
     * @return Const reference to the associated data
     */
    const T& operator[](Polytope::Type type) const
    {
      return *ptr(static_cast<size_t>(type));
    }

    /**
     * @brief Gets the number of polytope types (always constant).
     * @return Number of polytope types supported
     */
    constexpr size_t size() const { return Count; }

    /**
     * @brief Serialization method for Boost.Serialization.
     * @tparam Archive Archive type for serialization
     * @param ar Archive object
     * @param version Serialization version (unused)
     */
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
      for (size_t i = 0; i < Count; ++i)
        ar & *ptr(i);
    }

  private:
    /// Storage array for all polytope type mappings
    std::array<Storage, Count> m_map;

    /**
     * @brief Gets a pointer to the data at the given index.
     * @param i Index into the storage array
     * @return Pointer to the data object
     */
    T* ptr(size_t i)
    {
      assert(i < Count);
      return reinterpret_cast<T*>(&m_map[i]);
    }

    /**
     * @brief Gets a const pointer to the data at the given index.
     * @param i Index into the storage array
     * @return Const pointer to the data object
     */
    const T* ptr(size_t i) const
    {
      assert(i < Count);
      return reinterpret_cast<const T*>(&m_map[i]);
    }
  };
}

#endif
