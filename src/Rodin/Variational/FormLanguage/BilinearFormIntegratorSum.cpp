/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "UnaryMinus.h"
#include "BilinearFormIntegratorSum.h"

namespace Rodin::Variational::FormLanguage
{
   BilinearFormIntegratorSum::BilinearFormIntegratorSum(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorBase& rhs)
   {
      switch (lhs.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
            m_bfiDomainList.emplace_back(lhs.copy());
            break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }
      switch (rhs.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
            m_bfiDomainList.emplace_back(rhs.copy());
            break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }
   }

   BilinearFormIntegratorSum::BilinearFormIntegratorSum(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorBase& rhs)
   {
      m_bfiDomainList.reserve(lhs.m_bfiDomainList.size() + 1);
      for (const auto& p : lhs.m_bfiDomainList)
         m_bfiDomainList.emplace_back(p->copy());
      switch (rhs.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
            m_bfiDomainList.emplace_back(rhs.copy());
            break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }

   }

   BilinearFormIntegratorSum::BilinearFormIntegratorSum(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorSum& rhs)
   {
      m_bfiDomainList.reserve(lhs.m_bfiDomainList.size() + rhs.m_bfiDomainList.size());
      for (const auto& p : lhs.m_bfiDomainList)
         m_bfiDomainList.emplace_back(p->copy());
      for (const auto& p : rhs.m_bfiDomainList)
         m_bfiDomainList.emplace_back(p->copy());
   }

   BilinearFormIntegratorSum::BilinearFormIntegratorSum(
         const BilinearFormIntegratorSum& other)
   {
      m_bfiDomainList.reserve(other.m_bfiDomainList.size());
      for (const auto& p : other.m_bfiDomainList)
         m_bfiDomainList.emplace_back(p->copy());
   }

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorBase& rhs)
   {
      return BilinearFormIntegratorSum(lhs, rhs);
   }

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorBase& rhs)
   {
      return BilinearFormIntegratorSum(lhs, rhs);
   }

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorSum& rhs)
   {
      return BilinearFormIntegratorSum(rhs, lhs);
   }

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorSum& rhs)
   {
      return BilinearFormIntegratorSum(lhs, rhs);
   }

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorBase& rhs)
   {
      return BilinearFormIntegratorSum(
            lhs, UnaryMinus<BilinearFormIntegratorBase>(rhs));
   }

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorBase& rhs)
   {
      return BilinearFormIntegratorSum(
            lhs, UnaryMinus<BilinearFormIntegratorBase>(rhs));
   }

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorSum& rhs)
   {
      return BilinearFormIntegratorSum(
            rhs, UnaryMinus<BilinearFormIntegratorBase>(lhs));
   }

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorSum& rhs)
   {
      return BilinearFormIntegratorSum(
            lhs, UnaryMinus<BilinearFormIntegratorSum>(rhs));
   }
}

