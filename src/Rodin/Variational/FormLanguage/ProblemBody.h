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

#include "Rodin/Utility/OptionalReference.h"

#include "ForwardDecls.h"

#include "Base.h"
#include "BoundaryConditionList.h"

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Represents the body of a variational problem.
    */
   class ProblemBody : public Base
   {
      public:
         ProblemBody(
               const BilinearFormIntegratorBase& bfi);

         ProblemBody(
               const BilinearFormIntegratorBase& bfi,
               const LinearFormIntegratorBase& lfi);

         ProblemBody(
               const BilinearFormIntegratorBase& bfi,
               const BoundaryConditionList& bcs);

         ProblemBody(
               const BilinearFormIntegratorBase& bfi,
               const LinearFormIntegratorBase& lfi,
               const BoundaryConditionList& bcs);

         ProblemBody(const ProblemBody& other);

         BilinearFormIntegratorBase& getBilinearFormIntegrator();

         Utility::OptionalReference<LinearFormIntegratorBase>
         getLinearFormIntegrator();

         BoundaryConditionList& getBoundaryConditionList();

         virtual ProblemBody* copy() const noexcept override
         {
            return new ProblemBody(*this);
         }

      private:
         std::unique_ptr<BilinearFormIntegratorBase> m_bfi;
         std::unique_ptr<LinearFormIntegratorBase> m_lfi;
         std::unique_ptr<BoundaryConditionList> m_bcs;
   };

   ProblemBody operator+(
         const BilinearFormIntegratorBase& bfi, const LinearFormIntegratorBase& lfi);

   ProblemBody operator-(
         const BilinearFormIntegratorBase& bfi, const LinearFormIntegratorBase& lfi);

   ProblemBody operator+(
         const BilinearFormIntegratorBase& bfi, const BoundaryConditionList& bcs);

   ProblemBody operator+(
         const ProblemBody& pb, const BoundaryConditionList& bcs);
}

#endif
