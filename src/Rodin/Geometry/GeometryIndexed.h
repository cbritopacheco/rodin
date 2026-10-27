/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_GEOMETRYINDEXED_H
#define RODIN_GEOMETRY_GEOMETRYINDEXED_H

/**
 * @file
 * @brief Container for mapping polytope types to associated data.
 */

#include <array>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <boost/serialization/access.hpp>

#include "Rodin/Alert/Exception.h"

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
   * @see <a href="class_rodin_1_1_geometry_1_1_polytope.html#a1d1cfd8ffb84e947f82999c682b666a7">Polytope::Type</a>
   */
  template <class T>
  class GeometryIndexed
  {
    /// Number of different polytope types
    static constexpr size_t Count = Polytope::Types.size();

    friend class boost::serialization::access;

    public:
      /// @brief Default constructor that initializes all entries with default values.
      GeometryIndexed()
        : m_map{}
      {}

      /**
       * @brief Constructs from an initializer list of polytope type-value pairs.
       * @param l Initializer list of (polytope type, value) pairs
       * @note The list must contain exactly one entry for each polytope type
       */
      GeometryIndexed(std::initializer_list<std::pair<Polytope::Type, T>> l)
        : GeometryIndexed(l, std::make_index_sequence<Count>{})
      {}

      /**
       * @brief Copy constructor.
       * @param other GeometryIndexed object to copy from
       */
      GeometryIndexed(const GeometryIndexed& other) = default;

      /**
       * @brief Move constructor.
       * @param other GeometryIndexed object to move from
       */
      GeometryIndexed(GeometryIndexed&& other) noexcept(
        std::is_nothrow_move_constructible<T>::value) = default;

      /// @brief Destroys all retained values, including during constructor unwinding.
      ~GeometryIndexed() = default;

      /**
       * @brief Copy assignment operator.
       * @param other GeometryIndexed object to copy from
       * @return Reference to this object
       */
      GeometryIndexed& operator=(const GeometryIndexed& other) noexcept(
        std::is_nothrow_copy_assignable<T>::value)
      {
        if (this != &other)
        {
          for (size_t i = 0; i < Count; ++i)
            *ptr(i) = *other.ptr(i);
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param other GeometryIndexed object to move from
       * @return Reference to this object
       */
      GeometryIndexed& operator=(GeometryIndexed&& other) noexcept(
        std::is_nothrow_move_assignable<T>::value)
      {
        if (this != &other)
        {
          for (size_t i = 0; i < Count; ++i)
            *ptr(i) = std::move(*other.ptr(i));
        }
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
      constexpr size_t size() const
      {
        return Count;
      }

      /**
       * @brief Serialization method for Boost.Serialization.
       * @tparam Archive Archive type for serialization
       * @param ar Archive object
       * @param version Serialization version (unused)
       */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        for (size_t i = 0; i < Count; ++i)
          ar&* ptr(i);
      }

    private:
      /**
       * @brief Constructs entries in storage order with automatic exception cleanup.
       * @param entries Complete geometry-value mapping.
       * @param indices Compile-time storage indices.
       */
      template <size_t... Indices>
      GeometryIndexed(std::initializer_list<std::pair<Polytope::Type, T>> entries,
        [[maybe_unused]] std::index_sequence<Indices...> indices)
        : m_map{getEntry(entries, Indices)...}
      {}

      /**
       * @brief Finds a unique initializer entry for one geometry.
       * @param entries Complete geometry-value mapping.
       * @param index Geometry's storage index.
       * @returns The supplied value; missing and repeated geometries raise an exception.
       */
      static const T& getEntry(
        std::initializer_list<std::pair<Polytope::Type, T>> entries, size_t index)
      {
        const T* value = nullptr;
        for (const auto& entry : entries)
        {
          if (static_cast<size_t>(entry.first) == index)
          {
            if (value)
            {
              Alert::Exception()
                << "Repeated GeometryIndexed initializer geometry." << Alert::Raise;
            }
            value = &entry.second;
          }
        }
        if (entries.size() != Count || !value)
        {
          Alert::Exception() << "Incomplete GeometryIndexed initializer mapping."
                             << Alert::Raise;
        }
        return *value;
      }

      /// Storage array for all polytope type mappings
      std::array<T, Count> m_map;

      /**
       * @brief Gets a pointer to the data at the given index.
       * @param i Index into the storage array
       * @return Pointer to the data object
       */
      T* ptr(size_t i)
      {
        assert(i < Count);
        return &m_map[i];
      }

      /**
       * @brief Gets a const pointer to the data at the given index.
       * @param i Index into the storage array
       * @return Const pointer to the data object
       */
      const T* ptr(size_t i) const
      {
        assert(i < Count);
        return &m_map[i];
      }
  };
}

#endif
