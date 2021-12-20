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

#include "ProblemBody.h"

namespace Rodin::Variational::FormLanguage
{
   ProblemBody::ProblemBody(const BilinearFormIntegratorBase& bfi)
      :  m_bfi(bfi.copy()),
         m_bcs(std::make_unique<BoundaryConditionList>())
   {}

   ProblemBody::ProblemBody(
         const BilinearFormIntegratorBase& bfi,
         const LinearFormIntegratorBase& lfi)
      :  m_bfi(bfi.copy()),
         m_lfi(lfi.copy()),
         m_bcs(std::make_unique<BoundaryConditionList>())
   {}

   ProblemBody::ProblemBody(
         const BilinearFormIntegratorBase& bfi,
         const BoundaryConditionList& bcs)
      :  m_bfi(bfi.copy()),
         m_bcs(bcs.copy())
   {}

   ProblemBody::ProblemBody(
         const BilinearFormIntegratorBase& bfi,
         const LinearFormIntegratorBase& lfi,
         const BoundaryConditionList& bcs)
      :  m_bfi(bfi.copy()),
         m_lfi(lfi.copy()),
         m_bcs(bcs.copy())
   {}

   ProblemBody::ProblemBody(const ProblemBody& other)
      :  m_bfi(other.m_bfi->copy()),
         m_bcs(other.m_bcs->copy())
   {
      if (other.m_lfi)
         m_lfi = std::unique_ptr<LinearFormIntegratorBase>(other.m_lfi->copy());
   }

   BilinearFormIntegratorBase& ProblemBody::getBilinearFormIntegrator()
   {
      assert(m_bfi);
      return *m_bfi;
   }

   Utility::OptionalReference<LinearFormIntegratorBase>
   ProblemBody::getLinearFormIntegrator()
   {
      if (m_lfi)
         return *m_lfi;
      else
         return {};
   }

   BoundaryConditionList& ProblemBody::getBoundaryConditionList()
   {
      assert(m_bcs);
      return *m_bcs;
   }

   ProblemBody operator+(
         const BilinearFormIntegratorBase& bfi, const LinearFormIntegratorBase& lfi)
   {
      return ProblemBody(bfi, lfi);
   }

   ProblemBody operator-(
         const BilinearFormIntegratorBase& bfi, const LinearFormIntegratorBase& lfi)
   {
      return ProblemBody(bfi, LinearFormIntegratorUnaryMinus(lfi));
   }

   ProblemBody operator+(
         const BilinearFormIntegratorBase& bfi, const BoundaryConditionList& bcs)
   {
      return ProblemBody(bfi, bcs);
   }

   ProblemBody operator+(
         const ProblemBody& pb, const BoundaryConditionList& bcs)
   {
      ProblemBody res(pb);
      res.getBoundaryConditionList() += bcs;
      return res;
   }
}
