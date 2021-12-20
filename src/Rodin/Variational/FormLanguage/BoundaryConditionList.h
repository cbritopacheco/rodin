/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BOUNDARYCONDITIONLIST_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BOUNDARYCONDITIONLIST_H

#include <memory>
#include <vector>

#include "Rodin/Variational/BoundaryCondition.h"

#include "ForwardDecls.h"

#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   class BoundaryConditionList : public Base
   {
      public:
         class Iterator
         {
            public:
               using iterator_category
                  = std::input_iterator_tag;
               using pointer
                  = std::vector<std::unique_ptr<BoundaryConditionBase>>::iterator;
               using value_type
                  = BoundaryConditionBase&;
               using reference
                  = BoundaryConditionBase&;

               Iterator(pointer it) : m_ptr(it) {}

               reference operator*() const { return *(*m_ptr); }

               pointer operator->() { return m_ptr; }

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
                  = std::vector<std::unique_ptr<BoundaryConditionBase>>::const_iterator;
               using value_type
                  = BoundaryConditionBase&;
               using const_reference
                  = BoundaryConditionBase&;

               ConstIterator(const_pointer it) : m_ptr(it) {}

               const_reference operator*() const { return *(*m_ptr); }

               const_pointer operator->() { return m_ptr; }

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

         BoundaryConditionList() = default;

         BoundaryConditionList(const BoundaryConditionBase& bc);

         BoundaryConditionList(const BoundaryConditionList& other);

         BoundaryConditionList& operator+=(const BoundaryConditionList& other);

         size_t size() const;

         Iterator begin();
         Iterator end();

         ConstIterator begin() const;
         ConstIterator end() const;

         friend
         BoundaryConditionList operator+(
               const BoundaryConditionList& a, const BoundaryConditionList& b);

         BoundaryConditionList* copy() const noexcept override
         {
            return new BoundaryConditionList(*this);
         }
      private:
         std::vector<std::unique_ptr<BoundaryConditionBase>> m_bcList;
   };

}

#endif

