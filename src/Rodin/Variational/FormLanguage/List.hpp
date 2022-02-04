/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_LIST_HPP
#define RODIN_VARIATIONAL_FORMLANGUAGE_LIST_HPP

#include "List.h"

namespace Rodin::Variational::FormLanguage
{
   template <class T>
   List<T>::List(const T& v)
   {
      m_list.emplace_back(v.copy());
   }

   template <class T>
   List<T>::List(const List& other)
   {
      m_list.reserve(other.m_list.size());
      for (const auto& v : other.m_list)
         m_list.emplace_back(v->copy());
   }

   template <class T>
   typename List<T>::Iterator List<T>::begin()
   {
      return m_list.begin();
   }

   template <class T>
   typename List<T>::Iterator List<T>::end()
   {
      return m_list.end();
   }

   template <class T>
   typename List<T>::ConstIterator List<T>::begin() const
   {
      return m_list.begin();
   }

   template <class T>
   typename List<T>::ConstIterator List<T>::end() const
   {
      return m_list.end();
   }

   template <class T>
   size_t List<T>::size() const
   {
      return m_list.size();
   }

   template <class T>
   T& List<T>::append(const T& v)
   {
      return *m_list.emplace_back(v.copy());
   }

   template <class T>
   List<T> operator+(
         const List<T>& a, const List<T>& b)
   {
      List<T> res;
      res.m_bcList.reserve(a.m_bcList.size() + b.m_bcList.size());
      for (const auto& bc : a.m_bcList)
         res.m_bcList.emplace_back(bc->copy());
      for (const auto& bc : b.m_bcList)
         res.m_bcList.emplace_back(bc->copy());
      return res;
   }
}

#endif

