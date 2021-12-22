/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_PROBLEMBODY_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_PROBLEMBODY_H

#include <memory>
#include <optional>

#include "Rodin/Variational/BoundaryCondition.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"

#include "Base.h"
#include "List.h"


namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Represents the body of a variational problem.
    */
   class ProblemBody : public Base
   {
      public:
         ProblemBody(const BilinearFormDomainIntegrator& bfi);

         ProblemBody(const ProblemBody& other) = default;

         List<BoundaryConditionBase>& getBoundaryConditionList();
         List<LinearFormIntegratorBase>& getLinearFormDomainIntegratorList();
         List<LinearFormIntegratorBase>& getLinearFormBoundaryIntegratorList();
         List<BilinearFormIntegratorBase>& getBilinearFormDomainIntegratorList();

         virtual ProblemBody* copy() const noexcept override
         {
            return new ProblemBody(*this);
         }

      private:
         List<BoundaryConditionBase> m_bcList;
         List<BilinearFormIntegratorBase> m_bfiDomainList;
         List<LinearFormIntegratorBase> m_lfiDomainList;
         List<LinearFormIntegratorBase> m_lfiBoundaryList;
   };

   ProblemBody operator+(
         const ProblemBody& pb, const BilinearFormDomainIntegrator& bfi);

   ProblemBody operator-(
         const ProblemBody& pb, const BilinearFormDomainIntegrator& lfi);

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormDomainIntegrator& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormDomainIntegrator& lfi);

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormBoundaryIntegrator& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormBoundaryIntegrator& lfi);

   ProblemBody operator+(
         const ProblemBody& pb, const List<BoundaryConditionBase>& bcs);
}

#endif
