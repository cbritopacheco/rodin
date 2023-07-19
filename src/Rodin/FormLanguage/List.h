/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>
#include <memory>

#include "Base.h"

#ifndef RODIN_FORMLANGUAGE_LIST_H
#define RODIN_FORMLANGUAGE_LIST_H

namespace Rodin::FormLanguage
{
  template <class T>
  class List : public Base
  {
    public:
      using reference = T&;
      using const_reference = const T&;

      class Iterator
      {
        using internal_iterator = typename std::vector<std::unique_ptr<T>>::iterator;
        public:
          explicit constexpr Iterator(internal_iterator it) : m_it(it) {}
          constexpr Iterator& operator++() { m_it++; return *this; }
          constexpr Iterator operator++(int) { Iterator r = *this; ++(*this); return r; }
          constexpr bool operator==(const Iterator& other) const { return m_it == other.m_it; }
          constexpr bool operator!=(const Iterator& other) const { return !(*this == other); }
          constexpr reference operator*() const { assert(*m_it); return static_cast<reference>(**m_it); }
        private:
          internal_iterator m_it;
      };

      class ConstIterator
      {
        using internal_const_iterator = typename std::vector<std::unique_ptr<T>>::const_iterator;
        public:
          explicit constexpr ConstIterator(internal_const_iterator it) : m_it(it) {}
          constexpr ConstIterator& operator++() { m_it++; return *this; }
          constexpr ConstIterator operator++(int) { ConstIterator r = *this; ++(*this); return r; }
          constexpr bool operator==(const ConstIterator& other) const { return m_it == other.m_it; }
          constexpr bool operator!=(const ConstIterator& other) const { return !(*this == other); }
          constexpr const_reference operator*() const { assert(*m_it); return static_cast<const_reference>(**m_it); }
        private:
          internal_const_iterator m_it;
      };

      constexpr List() = default;

      constexpr List(const List& other)
        : Base(other)
      {
        m_list.reserve(other.m_list.size());
        for (const auto& p : other.m_list)
          m_list.emplace_back(p->copy());
      }

      constexpr List(List&& other)
        : Base(std::move(other)),
          m_list(std::move(other.m_list))
      {}

      virtual ~List() = default;

      constexpr List& operator=(const List& other)
      {
        if (this != &other)
        {
          m_list.reserve(other.m_list.size());
          for (const auto& p : other.m_list)
            m_list.emplace_back(p->copy());
        }
        return *this;
      }

      constexpr List& operator=(List&& other)
      {
        m_list = std::move(other.m_list);
        return *this;
      }

      constexpr reference at(size_t i)
      {
        return *m_list.at(i);
      }

      constexpr const_reference at(size_t i) const
      {
        return *m_list.at(i);
      }

      constexpr List& add(const T& v)
      {
        m_list.emplace_back(v.copy());
        return *this;
      }

      constexpr List& add(const FormLanguage::List<T>& other)
      {
        m_list.reserve(m_list.size() + other.m_list.size());
        for (const auto& p : other.m_list)
          m_list.emplace_back(p->copy());
        return *this;
      }

      constexpr List& clear()
      {
        m_list.clear();
        return *this;
      }

      constexpr size_t size() const
      {
        return m_list.size();
      }

      constexpr Iterator begin() noexcept
      {
        return Iterator(m_list.begin());
      }

      constexpr Iterator end() noexcept
      {
        return Iterator(m_list.end());
      }

      constexpr ConstIterator begin() const noexcept
      {
        return ConstIterator(m_list.begin());
      }

      constexpr ConstIterator end() const noexcept
      {
        return ConstIterator(m_list.end());
      }

      constexpr ConstIterator cbegin() const noexcept
      {
        return ConstIterator(m_list.cbegin());
      }

      constexpr ConstIterator cend() const noexcept
      {
        return ConstIterator(m_list.cend());
      }

      virtual List* copy() const noexcept override
      {
        return new List(*this);
      }

    private:
      std::vector<std::unique_ptr<T>> m_list;
  };
}

#endif
