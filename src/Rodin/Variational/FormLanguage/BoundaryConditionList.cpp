/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "BoundaryConditionList.h"

namespace Rodin::Variational::FormLanguage
{
   BoundaryConditionList::BoundaryConditionList(const BoundaryConditionBase& bc)
   {
      m_bcList.emplace_back(bc.copy());
   }

   BoundaryConditionList::BoundaryConditionList(const BoundaryConditionList& other)
   {
      m_bcList.reserve(other.m_bcList.size());
      for (const auto& bc : other.m_bcList)
         m_bcList.emplace_back(bc->copy());
   }

   BoundaryConditionList::Iterator BoundaryConditionList::begin()
   {
      return m_bcList.begin();
   }

   BoundaryConditionList::Iterator BoundaryConditionList::end()
   {
      return m_bcList.end();
   }

   BoundaryConditionList::ConstIterator BoundaryConditionList::begin() const
   {
      return m_bcList.begin();
   }

   BoundaryConditionList::ConstIterator BoundaryConditionList::end() const
   {
      return m_bcList.end();
   }

   size_t BoundaryConditionList::size() const
   {
      return m_bcList.size();
   }

   BoundaryConditionList& BoundaryConditionList::operator+=(
         const BoundaryConditionList& other)
   {
      m_bcList.reserve(m_bcList.size() + other.m_bcList.size());
      for (const auto& bc : other.m_bcList)
         m_bcList.emplace_back(bc->copy());
      return *this;
   }

   BoundaryConditionList operator+(
         const BoundaryConditionList& a, const BoundaryConditionList& b)
   {
      BoundaryConditionList res;
      res.m_bcList.reserve(a.m_bcList.size() + b.m_bcList.size());
      for (const auto& bc : a.m_bcList)
         res.m_bcList.emplace_back(bc->copy());
      for (const auto& bc : b.m_bcList)
         res.m_bcList.emplace_back(bc->copy());
      return res;
   }
}
