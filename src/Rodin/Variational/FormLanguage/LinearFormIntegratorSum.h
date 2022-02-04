/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATORSUM_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATORSUM_H

#include <memory>
#include <utility>

#include "Rodin/Variational/LinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   namespace Internal
   {
      class LinearFormIntegratorSum : public mfem::LinearFormIntegrator
      {
         public:
            LinearFormIntegratorSum(
                  mfem::LinearFormIntegrator& a, mfem::LinearFormIntegrator& b)
               : m_a(a), m_b(b)
            {}

            void AssembleRHSElementVect(
                  const mfem::FiniteElement& el,
                  mfem::ElementTransformation& Tr,
                  mfem::Vector& elvect)
            {
               m_a.AssembleRHSElementVect(el, Tr, elvect);
               mfem::Vector tmp;
               m_b.AssembleRHSElementVect(el, Tr, tmp);
               elvect.Add(1.0, tmp);
            }

         private:
            mfem::LinearFormIntegrator& m_a;
            mfem::LinearFormIntegrator& m_b;
      };
   }

   class LinearFormIntegratorSum : public LinearFormIntegratorBase
   {
      public:
         LinearFormIntegratorSum(
               const LinearFormIntegratorBase& lhs,
               const LinearFormIntegratorBase& rhs);

         LinearFormIntegratorSum(const LinearFormIntegratorSum& other);

         LinearFormIntegratorBase& getLHS();

         const LinearFormIntegratorBase& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         LinearFormIntegratorBase& getRHS();

         const LinearFormIntegratorBase& getRHS() const
         {
            assert(m_rhs);
            return *m_rhs;
         }

         const std::set<int>& getAttributes() const override
         {
            return getLHS().getAttributes();
         }

         void buildMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() override;

         LinearFormIntegratorSum* copy() const noexcept override
         {
            return new LinearFormIntegratorSum(*this);
         }

      private:
         std::unique_ptr<LinearFormIntegratorBase> m_lhs;
         std::unique_ptr<LinearFormIntegratorBase> m_rhs;
         std::unique_ptr<Internal::LinearFormIntegratorSum> m_mfemLFI;
   };

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs);

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs);
}

#endif

