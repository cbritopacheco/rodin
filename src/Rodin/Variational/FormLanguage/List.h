/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_LIST_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_LIST_H

#include <memory>
#include <vector>

#include "ForwardDecls.h"

#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   template <class T>
   class List : public Base
   {
      static_assert(std::is_base_of_v<Base, T>,
            "T must be derived from FormLanguage::Base");
      public:
         class Iterator
         {
            public:
               using iterator_category
                  = std::input_iterator_tag;
               using pointer
                  = typename std::vector<std::unique_ptr<T>>::iterator;
               using value_type = T&;
               using reference = T&;

               Iterator(pointer it) : m_ptr(it) {}

               reference operator*() const { return *(*m_ptr); }

               Iterator& operator++() { m_ptr++; return *this; }

               Iterator operator++(int)
               {
                  Iterator tmp = *this;
                  ++(*this);
                  return tmp;
               }

               friend
               bool operator==(const Iterator& a, const Iterator& b)
               {
                  return a.m_ptr == b.m_ptr;
               };

               friend
               bool operator!=(const Iterator& a, const Iterator& b)
               {
                  return a.m_ptr != b.m_ptr;
               };

            private:
               pointer m_ptr;
         };

         class ConstIterator
         {
            public:
               using iterator_category
                  = std::input_iterator_tag;
               using const_pointer
                  = typename std::vector<std::unique_ptr<T>>::const_iterator;
               using value_type = T&;
               using const_reference = const T&;

               ConstIterator(const_pointer it) : m_ptr(it) {}

               const_reference operator*() const { return *(*m_ptr); }

               ConstIterator& operator++() { m_ptr++; return *this; }

               ConstIterator operator++(int)
               {
                  ConstIterator tmp = *this;
                  ++(*this);
                  return tmp;
               }

               friend
               bool operator==(const ConstIterator& a, const ConstIterator& b)
               {
                  return a.m_ptr == b.m_ptr;
               };

               friend
               bool operator!=(const ConstIterator& a, const ConstIterator& b)
               {
                  return a.m_ptr != b.m_ptr;
               };

            private:
               const_pointer m_ptr;
         };

         List() = default;

         virtual ~List() = default;

         List(const T& bc);

         List(const List& other);

         List(List&& other) = default;

         List& operator=(const List& other)
         {
            if (this != &other)
            {
               List tmp(other);
               m_list.swap(tmp.m_list);
            }
            return *this;
         }

         List& operator+=(const List& other);

         T& append(const T& v);

         size_t size() const;

         Iterator begin();
         Iterator end();

         T& front()
         {
            return *(*m_list.begin());
         }

         T& back()
         {
            assert(size() >= 1);
            return *(*(m_list.end() - 1));
         }

         void clear() noexcept
         {
            m_list.clear();
         }

         ConstIterator begin() const;
         ConstIterator end() const;

         friend
         List operator+(const List& a, const List& b);

         List* copy() const noexcept override
         {
            return new List(*this);
         }
      private:
         std::vector<std::unique_ptr<T>> m_list;
   };
}

#include "List.hpp"

#endif
