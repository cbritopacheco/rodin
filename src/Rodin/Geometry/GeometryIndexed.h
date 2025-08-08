/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GEOMETRYINDEXED_H
#define RODIN_VARIATIONAL_GEOMETRYINDEXED_H

#include <array>
#include <cassert>
#include <cstddef>
#include <new>
#include <type_traits>
#include <boost/serialization/access.hpp>

#include "Polytope.h"

namespace Rodin::Geometry
{
  template <class T>
  class GeometryIndexed
  {
    static constexpr size_t Count = Polytope::Types.size();
    using Storage = typename std::aligned_storage<sizeof(T), alignof(T)>::type;

    friend class boost::serialization::access;

  public:
    GeometryIndexed()
    {
      for (size_t i = 0; i < Count; ++i)
        new (&m_map[i]) T();
    }

    GeometryIndexed(std::initializer_list<std::pair<Polytope::Type, T>> l)
    {
      assert(l.size() == Count);
      for (const auto& [type, value] : l)
        new (&m_map[static_cast<size_t>(type)]) T(std::move(value));
    }

    GeometryIndexed(const GeometryIndexed& other)
    {
      for (size_t i = 0; i < Count; ++i)
        new (&m_map[i]) T(*other.ptr(i));
    }

    GeometryIndexed(GeometryIndexed&& other) noexcept(std::is_nothrow_move_constructible<T>::value)
    {
      for (size_t i = 0; i < Count; ++i)
        new (&m_map[i]) T(std::move(*other.ptr(i)));
    }

    ~GeometryIndexed()
    {
      for (size_t i = 0; i < Count; ++i)
        ptr(i)->~T();
    }

    GeometryIndexed& operator=(const GeometryIndexed& other)
      noexcept(std::is_nothrow_copy_assignable<T>::value)
    {
      if (this != &other)
        for (size_t i = 0; i < Count; ++i)
          *ptr(i) = *other.ptr(i);
      return *this;
    }

    GeometryIndexed& operator=(GeometryIndexed&& other)
      noexcept(std::is_nothrow_move_assignable<T>::value)
    {
      if (this != &other)
        for (size_t i = 0; i < Count; ++i)
          *ptr(i) = std::move(*other.ptr(i));
      return *this;
    }

    T& operator[](Polytope::Type type)
    {
      return *ptr(static_cast<size_t>(type));
    }

    const T& operator[](Polytope::Type type) const
    {
      return *ptr(static_cast<size_t>(type));
    }

    constexpr size_t size() const { return Count; }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
      for (size_t i = 0; i < Count; ++i)
        ar & *ptr(i);
    }

  private:
    std::array<Storage, Count> m_map;

    T* ptr(size_t i)
    {
      assert(i < Count);
      return reinterpret_cast<T*>(&m_map[i]);
    }

    const T* ptr(size_t i) const
    {
      assert(i < Count);
      return reinterpret_cast<const T*>(&m_map[i]);
    }
  };
}

#endif
