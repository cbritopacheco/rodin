/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file List.h
 * @brief Polymorphic container for form language objects.
 *
 * This file defines the List class template, which provides a type-safe
 * container for storing polymorphic form language objects. The List maintains
 * ownership of its elements through unique pointers and supports standard
 * container operations.
 *
 * ## Design
 * The List class is designed to work with types derived from FormLanguage::Base,
 * leveraging the copy() method for proper polymorphic copying. Elements are
 * stored as unique_ptr instances, ensuring automatic memory management.
 *
 * ## Usage Example
 * @code{.cpp}
 * List<BoundaryCondition> bcList;
 * bcList.add(DirichletBC(...));
 * bcList.add(NeumannBC(...));
 * 
 * for (const auto& bc : bcList)
 * {
 *   // Process each boundary condition
 * }
 * @endcode
 */
#include <vector>
#include <memory>

#include "Base.h"

#ifndef RODIN_FORMLANGUAGE_LIST_H
#define RODIN_FORMLANGUAGE_LIST_H

namespace Rodin::FormLanguage
{
  /**
   * @brief Polymorphic container for form language objects.
   * @tparam T Element type (must be derived from or compatible with Base)
   *
   * List provides a dynamic container for storing polymorphic form language
   * objects with value semantics. The container manages object lifetimes
   * through unique pointers and provides deep copying capabilities.
   *
   * ## Features
   * - **Polymorphic Storage**: Stores derived types through base class pointers
   * - **Value Semantics**: Copy constructor and assignment perform deep copies
   * - **Iterator Support**: Provides standard forward iterators
   * - **Type Safety**: Elements are accessed through their stored type T
   *
   * ## Memory Management
   * - Elements are stored as std::unique_ptr<T>
   * - Copy operations use the polymorphic copy() method
   * - Automatic cleanup on destruction
   *
   * ## Complexity
   * - Access: @f$ O(1) @f$ for at()
   * - Insertion: @f$ O(1) @f$ amortized for add()
   * - Copy: @f$ O(n) @f$ where @f$ n @f$ is the number of elements
   */
  template <class T>
  class List : public Base
  {
    public:
      /**
       * @brief Reference type for list elements.
       */
      using reference = T&;

      /**
       * @brief Const reference type for list elements.
       */
      using const_reference = const T&;

      /**
       * @brief Forward iterator for non-const List.
       *
       * Iterator provides forward iteration over the list elements,
       * automatically dereferencing the unique_ptr storage to access
       * the contained objects.
       */
      class Iterator
      {
        using internal_iterator = typename std::vector<std::unique_ptr<T>>::iterator;
        public:
          /**
           * @brief Constructs an iterator from an internal iterator.
           * @param[in] it Internal vector iterator
           */
          explicit constexpr Iterator(internal_iterator it) : m_it(it) {}

          /**
           * @brief Pre-increment operator.
           * @return Reference to this iterator after increment
           */
          constexpr Iterator& operator++() { m_it++; return *this; }

          /**
           * @brief Post-increment operator.
           * @return Copy of iterator before increment
           */
          constexpr Iterator operator++(int) { Iterator r = *this; ++(*this); return r; }

          /**
           * @brief Equality comparison operator.
           * @param[in] other Iterator to compare with
           * @return True if iterators point to the same element
           */
          constexpr bool operator==(const Iterator& other) const { return m_it == other.m_it; }

          /**
           * @brief Inequality comparison operator.
           * @param[in] other Iterator to compare with
           * @return True if iterators point to different elements
           */
          constexpr bool operator!=(const Iterator& other) const { return !(*this == other); }

          /**
           * @brief Dereference operator.
           * @return Reference to the element at the current position
           */
          constexpr reference operator*() const { assert(*m_it); return static_cast<reference>(**m_it); }
        private:
          internal_iterator m_it;
      };

      /**
       * @brief Forward iterator for const List.
       *
       * ConstIterator provides forward iteration over const list elements,
       * ensuring read-only access to the contained objects.
       */
      class ConstIterator
      {
        using internal_const_iterator = typename std::vector<std::unique_ptr<T>>::const_iterator;
        public:
          /**
           * @brief Constructs a const iterator from an internal const iterator.
           * @param[in] it Internal vector const iterator
           */
          explicit constexpr ConstIterator(internal_const_iterator it) : m_it(it) {}

          /**
           * @brief Pre-increment operator.
           * @return Reference to this iterator after increment
           */
          constexpr ConstIterator& operator++() { m_it++; return *this; }

          /**
           * @brief Post-increment operator.
           * @return Copy of iterator before increment
           */
          constexpr ConstIterator operator++(int) { ConstIterator r = *this; ++(*this); return r; }

          /**
           * @brief Equality comparison operator.
           * @param[in] other Iterator to compare with
           * @return True if iterators point to the same element
           */
          constexpr bool operator==(const ConstIterator& other) const { return m_it == other.m_it; }

          /**
           * @brief Inequality comparison operator.
           * @param[in] other Iterator to compare with
           * @return True if iterators point to different elements
           */
          constexpr bool operator!=(const ConstIterator& other) const { return !(*this == other); }

          /**
           * @brief Dereference operator.
           * @return Const reference to the element at the current position
           */
          constexpr const_reference operator*() const { assert(*m_it); return static_cast<const_reference>(**m_it); }
        private:
          internal_const_iterator m_it;
      };

      /**
       * @brief Default constructor creates an empty list.
       */
      constexpr List() = default;

      /**
       * @brief Copy constructor performs deep copy of all elements.
       * @param[in] other List to copy from
       *
       * Each element is copied using its polymorphic copy() method,
       * ensuring proper copying of derived types.
       */
      constexpr List(const List& other)
        : Base(other)
      {
        m_list.reserve(other.m_list.size());
        for (const auto& p : other.m_list)
          m_list.emplace_back(p->copy());
      }

      /**
       * @brief Move constructor transfers ownership of elements.
       * @param[in] other List to move from
       */
      constexpr List(List&& other)
        : Base(std::move(other)),
          m_list(std::move(other.m_list))
      {}

      /**
       * @brief Destructor.
       */
      virtual ~List() = default;

      /**
       * @brief Copy assignment operator performs deep copy.
       * @param[in] other List to copy from
       * @return Reference to this list
       *
       * Clears existing elements and performs a deep copy of all
       * elements from other using polymorphic copy().
       */
      constexpr List& operator=(const List& other)
      {
        if (this != &other)
        {
          m_list.clear();
          m_list.reserve(other.m_list.size());
          for (const auto& p : other.m_list)
            m_list.emplace_back(p->copy());
        }
        return *this;
      }

      /**
       * @brief Move assignment operator transfers ownership.
       * @param[in] other List to move from
       * @return Reference to this list
       */
      constexpr List& operator=(List&& other)
      {
        m_list = std::move(other.m_list);
        return *this;
      }

      /**
       * @brief Accesses element at specified index with bounds checking.
       * @param[in] i Index of element to access
       * @return Reference to element at index i
       * @throws std::out_of_range if i >= size()
       */
      constexpr reference at(size_t i)
      {
        return *m_list.at(i);
      }

      /**
       * @brief Accesses element at specified index with bounds checking (const version).
       * @param[in] i Index of element to access
       * @return Const reference to element at index i
       * @throws std::out_of_range if i >= size()
       */
      constexpr const_reference at(size_t i) const
      {
        return *m_list.at(i);
      }

      /**
       * @brief Adds a single element to the list.
       * @param[in] v Element to add (will be copied)
       * @return Reference to this list for chaining
       *
       * Creates a polymorphic copy of the element using copy()
       * and appends it to the list.
       */
      constexpr List& add(const T& v)
      {
        m_list.emplace_back(v.copy());
        return *this;
      }

      /**
       * @brief Adds all elements from another list.
       * @param[in] other List containing elements to add
       * @return Reference to this list for chaining
       *
       * Performs a deep copy of all elements from other and
       * appends them to this list.
       */
      constexpr List& add(const FormLanguage::List<T>& other)
      {
        m_list.reserve(m_list.size() + other.m_list.size());
        for (const auto& p : other.m_list)
          m_list.emplace_back(p->copy());
        return *this;
      }

      /**
       * @brief Removes all elements from the list.
       * @return Reference to this list for chaining
       */
      constexpr List& clear()
      {
        m_list.clear();
        return *this;
      }

      /**
       * @brief Returns the number of elements in the list.
       * @return Number of elements
       */
      constexpr size_t size() const
      {
        return m_list.size();
      }

      /**
       * @brief Checks if the list is empty.
       * @return True if the list contains no elements
       */
      constexpr bool empty() const
      {
        return m_list.empty();
      }

      /**
       * @brief Returns an iterator to the beginning.
       * @return Iterator to the first element
       */
      constexpr Iterator begin() noexcept
      {
        return Iterator(m_list.begin());
      }

      /**
       * @brief Returns an iterator to the end.
       * @return Iterator to the element following the last element
       */
      constexpr Iterator end() noexcept
      {
        return Iterator(m_list.end());
      }

      /**
       * @brief Returns a const iterator to the beginning.
       * @return Const iterator to the first element
       */
      constexpr ConstIterator begin() const noexcept
      {
        return ConstIterator(m_list.begin());
      }

      /**
       * @brief Returns a const iterator to the end.
       * @return Const iterator to the element following the last element
       */
      constexpr ConstIterator end() const noexcept
      {
        return ConstIterator(m_list.end());
      }

      /**
       * @brief Returns a const iterator to the beginning.
       * @return Const iterator to the first element
       */
      constexpr ConstIterator cbegin() const noexcept
      {
        return ConstIterator(m_list.cbegin());
      }

      /**
       * @brief Returns a const iterator to the end.
       * @return Const iterator to the element following the last element
       */
      constexpr ConstIterator cend() const noexcept
      {
        return ConstIterator(m_list.cend());
      }

      /**
       * @brief Creates a polymorphic copy of this list.
       * @return Non-owning pointer to the copied list
       *
       * Implements the Copyable interface, creating a deep copy
       * of the list and all its elements.
       */
      virtual List* copy() const noexcept override
      {
        return new List(*this);
      }

    private:
      /**
       * @brief Internal storage for list elements.
       */
      std::vector<std::unique_ptr<T>> m_list;
  };
}

#endif
