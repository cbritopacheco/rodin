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
#include "Rodin/Variational/EssentialBoundary.h"

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
         using LFIList = std::vector<std::unique_ptr<LinearFormIntegratorBase>>;
         using BFIList = std::vector<std::unique_ptr<BilinearFormIntegratorBase>>;

         ProblemBody(const BilinearFormIntegratorBase& bfi);

         ProblemBody(const BilinearFormIntegratorSum& bfi);

         ProblemBody(const ProblemBody& other);

         EssentialBoundary& getEssentialBoundary();

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
         EssentialBoundary m_essBdr;
   };

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

   ProblemBody operator+(
         const ProblemBody& pb, const LinearFormIntegratorSum& lfi);

   ProblemBody operator-(
         const ProblemBody& pb, const LinearFormIntegratorSum& lfi);

   template <class T, class Value>
   ProblemBody operator+(const ProblemBody& pb, const DirichletBC<T, Value>& bc)
   {
      ProblemBody res(pb);
      res.getEssentialBoundary().add(bc);
      return res;
   }
}

#endif
