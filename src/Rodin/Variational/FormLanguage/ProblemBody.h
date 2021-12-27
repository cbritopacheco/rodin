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

#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/NeumannBC.h"

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
         ProblemBody(const BilinearFormIntegratorBase& bfi);

         ProblemBody(const ProblemBody& other) = default;

         List<BoundaryConditionBase>& getDirichletBCList();

         List<BoundaryConditionBase>& getNeumannBCList();

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

         List<BoundaryConditionBase> m_dbcs;

         List<BoundaryConditionBase> m_nbcs;
   };

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormDomainIntegrator& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormDomainIntegrator& lfi);

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormBoundaryIntegrator& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormBoundaryIntegrator& lfi);

   template <class T>
   ProblemBody operator+(
         const ProblemBody& pb, const DirichletBC<T>& bc)
   {
      ProblemBody res(pb);
      res.getDirichletBCList().append(bc);
      return res;
   }

   template <class T>
   ProblemBody operator+(
         const ProblemBody& pb, const NeumannBC<T>& bc)
   {
      ProblemBody res(pb);
      res.getNeumannBCList().append(bc);
      return res;
   }
}

#endif
