/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEMBODY_H
#define RODIN_VARIATIONAL_PROBLEMBODY_H

#include <vector>
#include <memory>
#include <optional>

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"

#include "BilinearFormIntegrator.h"
#include "LinearFormIntegrator.h"
#include "EssentialBoundary.h"


namespace Rodin::Variational
{
   /**
    * @brief Represents the body of a variational problem.
    */
   class ProblemBody : public FormLanguage::Base
   {
      public:
         ProblemBody() = default;

         ProblemBody(const BilinearFormIntegratorBase& bfi)
         {
            m_bfis.add(bfi);
         }

         ProblemBody(const FormLanguage::List<BilinearFormIntegratorBase>& bfis)
            : m_bfis(bfis)
         {}

         ProblemBody(const ProblemBody& other)
            :  FormLanguage::Base(other),
               m_bfis(other.m_bfis),
               m_lfis(other.m_lfis),
               m_essBdr(other.m_essBdr)
         {}

         ProblemBody(ProblemBody&& other)
            :  FormLanguage::Base(std::move(other)),
               m_bfis(std::move(other.m_bfis)),
               m_lfis(std::move(other.m_lfis)),
               m_essBdr(std::move(other.m_essBdr))
         {}

         ProblemBody& operator=(ProblemBody&& other)
         {
            m_bfis = std::move(other.m_bfis);
            m_lfis = std::move(other.m_lfis);
            m_essBdr = std::move(other.m_essBdr);
            return *this;
         }

         EssentialBoundary& getEssentialBoundary()
         {
            return m_essBdr;
         }

         const EssentialBoundary& getEssentialBoundary() const
         {
            return m_essBdr;
         }

         FormLanguage::List<BilinearFormIntegratorBase>& getBFIs()
         {
            return m_bfis;
         }

         const FormLanguage::List<BilinearFormIntegratorBase>& getBFIs() const
         {
            return m_bfis;
         }

         FormLanguage::List<LinearFormIntegratorBase>& getLFIs()
         {
            return m_lfis;
         }

         const FormLanguage::List<LinearFormIntegratorBase>& getLFIs() const
         {
            return m_lfis;
         }

         virtual ProblemBody* copy() const noexcept override
         {
            return new ProblemBody(*this);
         }

      private:
         FormLanguage::List<BilinearFormIntegratorBase> m_bfis;
         FormLanguage::List<LinearFormIntegratorBase>   m_lfis;
         EssentialBoundary m_essBdr;
   };

   ProblemBody operator+(
      const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

   ProblemBody operator-(
      const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

   ProblemBody operator+(
      const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfi);

   ProblemBody operator-(
      const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfi);

   template <class T>
   ProblemBody operator+(const ProblemBody& pb, const DirichletBC<T>& bc)
   {
      ProblemBody res(pb);
      res.getEssentialBoundary().add(bc);
      return res;
   }
}

#endif
