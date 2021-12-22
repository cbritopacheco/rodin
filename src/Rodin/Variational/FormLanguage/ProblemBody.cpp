/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BoundaryCondition.h"

#include "LinearFormIntegratorUnaryMinus.h"
#include "BilinearFormIntegratorUnaryMinus.h"

#include "ProblemBody.h"

namespace Rodin::Variational::FormLanguage
{
   ProblemBody::ProblemBody(const BilinearFormDomainIntegrator& bfi)
      : m_bfiDomainList(bfi)
   {}

   List<BoundaryConditionBase>&
   ProblemBody::getBoundaryConditionList()
   {
      return m_bcList;
   }

   List<LinearFormIntegratorBase>&
   ProblemBody::getLinearFormDomainIntegratorList()
   {
      return m_lfiDomainList;
   }

   List<LinearFormIntegratorBase>&
   ProblemBody::getLinearFormBoundaryIntegratorList()
   {
      return m_lfiBoundaryList;
   }

   List<BilinearFormIntegratorBase>&
   ProblemBody::getBilinearFormDomainIntegratorList()
   {
      return m_bfiDomainList;
   }

   ProblemBody operator+(const ProblemBody& pb, const List<BoundaryConditionBase>& bcs)
   {
      ProblemBody res(pb);
      res.getBoundaryConditionList() += bcs;
      return res;
   }

   ProblemBody operator+(
         const ProblemBody& pb, const BilinearFormDomainIntegrator& bfi)
   {
      ProblemBody res(pb);
      res.getBilinearFormDomainIntegratorList() += bfi;
      return res;
   }

   ProblemBody operator-(
         const ProblemBody& pb, const BilinearFormDomainIntegrator& bfi)
   {
      ProblemBody res(pb);
      res.getBilinearFormDomainIntegratorList() += -bfi;
      return res;
   }

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormDomainIntegrator& lfi)
   {
      ProblemBody res(pb);
      // Sign is opposite because we want the LinearFormIntegrator on the LHS
      res.getLinearFormDomainIntegratorList() += -lfi;
      return res;
   }

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormDomainIntegrator& lfi)
   {
      ProblemBody res(pb);
      // Sign is opposite because we want the LinearFormIntegrator on the LHS
      res.getLinearFormDomainIntegratorList() += lfi;
      return res;
   }

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormBoundaryIntegrator& lfi)
   {
      ProblemBody res(pb);
      // Sign is opposite because we want the LinearFormIntegrator on the LHS
      res.getLinearFormBoundaryIntegratorList() += -lfi;
      return res;
   }

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormBoundaryIntegrator& lfi)
   {
      ProblemBody res(pb);
      // Sign is opposite because we want the LinearFormIntegrator on the LHS
      res.getLinearFormBoundaryIntegratorList() += lfi;
      return res;
   }
}
