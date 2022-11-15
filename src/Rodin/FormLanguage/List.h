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
      static_assert(std::is_base_of_v<Base, T>);

      public:
         using reference = T&;
         using const_reference = const T&;

         class iterator
         {
            using internal_iterator = typename std::vector<std::unique_ptr<T>>::iterator;

            public:
               explicit constexpr iterator(internal_iterator it) : m_it(it) {}
               constexpr iterator& operator++()
               { m_it++; return *this; }
               constexpr iterator operator++(int)
               { iterator r = *this; ++(*this); return r; }
               constexpr bool operator==(iterator other) const
               { return m_it == other.m_it; }
               constexpr bool operator!=(iterator other) const
               { return !(*this == other); }
               constexpr reference operator*() const
               { assert(*m_it); return static_cast<reference>(**m_it); }

            private:
               internal_iterator m_it;
         };

         class const_iterator
         {
            using internal_const_iterator = typename std::vector<std::unique_ptr<T>>::const_iterator;

            public:
               explicit constexpr const_iterator(internal_const_iterator it) : m_it(it) {}
               constexpr const_iterator& operator++()
               { m_it++; return *this; }
               constexpr const_iterator operator++(int)
               { const_iterator r = *this; ++(*this); return r; }
               constexpr bool operator==(const_iterator other) const
               { return m_it == other.m_it; }
               constexpr bool operator!=(const_iterator other) const
               { return !(*this == other); }
               constexpr const_reference operator*() const
               { assert(*m_it); return static_cast<const_reference>(**m_it); }

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

         constexpr iterator begin() noexcept
         {
            return iterator(m_list.begin());
         }

         constexpr iterator end() noexcept
         {
            return iterator(m_list.end());
         }

         constexpr const_iterator begin() const noexcept
         {
            return const_iterator(m_list.begin());
         }

         constexpr const_iterator end() const noexcept
         {
            return const_iterator(m_list.end());
         }

         constexpr const_iterator cbegin() const noexcept
         {
            return const_iterator(m_list.cbegin());
         }

         constexpr const_iterator cend() const noexcept
         {
            return const_iterator(m_list.cend());
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
