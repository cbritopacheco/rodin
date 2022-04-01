/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "../UnaryMinus.h"

#include "LinearFormIntegratorSum.h"

namespace Rodin::Variational::FormLanguage
{
   LinearFormIntegratorSum::LinearFormIntegratorSum(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs)
   {
      switch (lhs.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
            m_lfiDomainList.emplace_back(lhs.copy());
            break;
         case IntegratorRegion::Boundary:
            m_lfiBoundaryList.emplace_back(lhs.copy());
            break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }
      switch (rhs.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
            m_lfiDomainList.emplace_back(rhs.copy());
            break;
         case IntegratorRegion::Boundary:
            m_lfiBoundaryList.emplace_back(rhs.copy());
            break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }
   }

   LinearFormIntegratorSum::LinearFormIntegratorSum(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorBase& rhs)
   {
      m_lfiDomainList.reserve(lhs.m_lfiDomainList.size() + 1);
      for (const auto& p : lhs.m_lfiDomainList)
         m_lfiDomainList.emplace_back(p->copy());
      m_lfiBoundaryList.reserve(lhs.m_lfiBoundaryList.size() + 1);
      for (const auto& p : lhs.m_lfiBoundaryList)
         m_lfiBoundaryList.emplace_back(p->copy());
      switch (rhs.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
            m_lfiDomainList.emplace_back(rhs.copy());
            break;
         case IntegratorRegion::Boundary:
            m_lfiBoundaryList.emplace_back(rhs.copy());
            break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }

   }

   LinearFormIntegratorSum::LinearFormIntegratorSum(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorSum& rhs)
   {
      m_lfiDomainList.reserve(lhs.m_lfiDomainList.size() + rhs.m_lfiDomainList.size());
      for (const auto& p : lhs.m_lfiDomainList)
         m_lfiDomainList.emplace_back(p->copy());
      for (const auto& p : rhs.m_lfiDomainList)
         m_lfiDomainList.emplace_back(p->copy());
      m_lfiBoundaryList.reserve(lhs.m_lfiBoundaryList.size() + rhs.m_lfiBoundaryList.size());
      for (const auto& p : lhs.m_lfiBoundaryList)
         m_lfiBoundaryList.emplace_back(p->copy());
      for (const auto& p : rhs.m_lfiBoundaryList)
         m_lfiBoundaryList.emplace_back(p->copy());
   }

   LinearFormIntegratorSum::LinearFormIntegratorSum(
         const LinearFormIntegratorSum& other)
   {
      m_lfiDomainList.reserve(other.m_lfiDomainList.size());
      for (const auto& p : other.m_lfiDomainList)
         m_lfiDomainList.emplace_back(p->copy());
      m_lfiBoundaryList.reserve(other.m_lfiBoundaryList.size());
      for (const auto& p : other.m_lfiBoundaryList)
         m_lfiBoundaryList.emplace_back(p->copy());
   }

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs)
   {
      return LinearFormIntegratorSum(lhs, rhs);
   }

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorBase& rhs)
   {
      return LinearFormIntegratorSum(lhs, rhs);
   }

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorSum& rhs)
   {
      return LinearFormIntegratorSum(rhs, lhs);
   }

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorSum& rhs)
   {
      return LinearFormIntegratorSum(lhs, rhs);
   }

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs)
   {
      return LinearFormIntegratorSum(
            lhs, UnaryMinus<LinearFormIntegratorBase>(rhs));
   }

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorBase& rhs)
   {
      return LinearFormIntegratorSum(
            lhs, UnaryMinus<LinearFormIntegratorBase>(rhs));
   }

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorSum& rhs)
   {
      return LinearFormIntegratorSum(
            rhs, UnaryMinus<LinearFormIntegratorBase>(lhs));
   }

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorSum& rhs)
   {
      return LinearFormIntegratorSum(
            lhs, UnaryMinus<LinearFormIntegratorSum>(rhs));
   }
}

