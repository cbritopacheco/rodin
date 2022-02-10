/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "LinearFormIntegratorUnaryMinus.h"
#include "BilinearFormIntegratorUnaryMinus.h"

#include "ProblemBody.h"

namespace Rodin::Variational::FormLanguage
{
   ProblemBody::ProblemBody(const ProblemBody& other)
   {
      m_nbcs.reserve(other.m_nbcs.size());
      m_dbcs.reserve(other.m_dbcs.size());
      m_bfiDomainList.reserve(other.m_bfiDomainList.size());
      m_lfiDomainList.reserve(other.m_lfiDomainList.size());
      m_lfiBoundaryList.reserve(other.m_lfiBoundaryList.size());

      for (const auto& v : other.m_nbcs)
         m_nbcs.emplace_back(v->copy());
      for (const auto& v : other.m_dbcs)
         m_dbcs.emplace_back(v->copy());
      for (const auto& v : other.m_bfiDomainList)
         m_bfiDomainList.emplace_back(v->copy());
      for (const auto& v : other.m_lfiDomainList)
         m_lfiDomainList.emplace_back(v->copy());
      for (const auto& v : other.m_lfiBoundaryList)
         m_lfiBoundaryList.emplace_back(v->copy());
   }

   ProblemBody::ProblemBody(const BilinearFormIntegratorBase& bfi)
   {
      m_bfiDomainList.emplace_back(bfi.copy());
   }

   ProblemBody::BCList& ProblemBody::getNeumannBCList()
   {
      return m_nbcs;
   }

   ProblemBody::BCList& ProblemBody::getDirichletBCList()
   {
      return m_dbcs;
   }

   ProblemBody::LFIList&
   ProblemBody::getLinearFormDomainIntegratorList()
   {
      return m_lfiDomainList;
   }

   ProblemBody::LFIList&
   ProblemBody::getLinearFormBoundaryIntegratorList()
   {
      return m_lfiBoundaryList;
   }

   ProblemBody::BFIList&
   ProblemBody::getBilinearFormDomainIntegratorList()
   {
      return m_bfiDomainList;
   }

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormIntegratorBase& lfi)
   {
      ProblemBody res(pb);
      // Sign is opposite because we want the LinearFormIntegrator on the LHS
      switch (lfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         res.getLinearFormDomainIntegratorList().emplace_back(
               new LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>(lfi));
         break;
         case IntegratorRegion::Boundary:
         res.getLinearFormBoundaryIntegratorList().emplace_back(
               new LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>(lfi));
         break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }
      return res;
   }

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormIntegratorBase& lfi)
   {
      ProblemBody res(pb);
      // Sign is opposite because we want the LinearFormIntegrator on the LHS
      switch (lfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         res.getLinearFormDomainIntegratorList().emplace_back(lfi.copy());
         break;
         case IntegratorRegion::Boundary:
         res.getLinearFormBoundaryIntegratorList().emplace_back(lfi.copy());
         break;
         default:
            Alert::Exception() << "IntegratorRegion not supported" << Alert::Raise;
      }
      return res;
   }

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormIntegratorSum& lfi)
   {
      ProblemBody res(pb);
      for (const auto& p : lfi.getLinearFormDomainIntegratorList())
         res.getLinearFormDomainIntegratorList().emplace_back(p->copy());
      for (const auto& p : lfi.getLinearFormBoundaryIntegratorList())
         res.getLinearFormBoundaryIntegratorList().emplace_back(p->copy());
      return res;
   }

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormIntegratorSum& lfi)
   {
      ProblemBody res(pb);
      for (const auto& p : lfi.getLinearFormDomainIntegratorList())
         res.getLinearFormDomainIntegratorList().emplace_back(p->copy());
      for (const auto& p : lfi.getLinearFormBoundaryIntegratorList())
         res.getLinearFormBoundaryIntegratorList().emplace_back(p->copy());
      return res;
   }
}
