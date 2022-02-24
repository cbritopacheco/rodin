/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_PROBLEMBODY_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_PROBLEMBODY_H

#include <vector>
#include <memory>
#include <optional>

#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/DirichletBC.h"

#include "ForwardDecls.h"

#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Represents the body of a variational problem.
    */
   class ProblemBody : public Base
   {
      public:
         using BCList =  std::vector<DirichletBC>;
         using LFIList = std::vector<std::unique_ptr<LinearFormIntegratorBase>>;
         using BFIList = std::vector<std::unique_ptr<BilinearFormIntegratorBase>>;

         ProblemBody(const BilinearFormIntegratorBase& bfi);

         ProblemBody(const BilinearFormIntegratorSum& bfi);

         ProblemBody(const ProblemBody& other);

         BCList& getDirichletBCList();
         LFIList& getLinearFormDomainIntegratorList();
         LFIList& getLinearFormBoundaryIntegratorList();
         BFIList& getBilinearFormDomainIntegratorList();

         virtual ProblemBody* copy() const noexcept override
         {
            return new ProblemBody(*this);
         }

      private:
         BFIList m_bfiDomainList;
         LFIList m_lfiDomainList;
         LFIList m_lfiBoundaryList;

         BCList m_dbcs;
         BCList m_nbcs;
   };

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormIntegratorSum& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormIntegratorSum& lfi);

   ProblemBody operator+(const ProblemBody& pb, const DirichletBC& bc);
}

#endif
