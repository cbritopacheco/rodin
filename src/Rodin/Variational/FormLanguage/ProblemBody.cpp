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
   ProblemBody::ProblemBody(const BilinearFormIntegratorBase& bfi)
      : m_bfiDomainList(bfi)
   {}

   List<BoundaryConditionBase>& ProblemBody::getNeumannBCList()
   {
      return m_nbcs;
   }

   List<BoundaryConditionBase>& ProblemBody::getDirichletBCList()
   {
      return m_dbcs;
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
