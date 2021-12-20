/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATORUNARYMINUS_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATORUNARYMINUS_H

#include "Rodin/Variational/LinearFormIntegrator.h"

namespace Rodin::Variational::FormLanguage
{
   namespace Internal
   {
      class LFIUnaryMinus : public mfem::LinearFormIntegrator
      {
         public:
            LFIUnaryMinus(mfem::LinearFormIntegrator& lfi)
               : m_lfi(lfi)
            {}

            void AssembleRHSElementVect(
                  const mfem::FiniteElement& el,
                  mfem::ElementTransformation& Tr,
                  mfem::Vector& elvect)
            {
               m_lfi.AssembleRHSElementVect(el, Tr, elvect);
               elvect.Neg();
            }

         private:
            mfem::LinearFormIntegrator& m_lfi;
      };
   }

   /**
    * @internal
    */
   class LinearFormIntegratorUnaryMinus : public LinearFormIntegratorBase
   {
      public:
         LinearFormIntegratorUnaryMinus(const LinearFormIntegratorBase& lfi);

         LinearFormIntegratorUnaryMinus(const LinearFormIntegratorUnaryMinus& other);

         LinearFormIntegratorBase& getLFI();

         void buildMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() override;

         LinearFormIntegratorUnaryMinus* copy() const noexcept override
         {
            return new LinearFormIntegratorUnaryMinus(*this);
         }

      private:
         std::unique_ptr<LinearFormIntegratorBase> m_lfi;
         std::unique_ptr<Internal::LFIUnaryMinus> m_mfemLFI;
   };

   LinearFormIntegratorUnaryMinus operator-(const LinearFormIntegratorBase& lfi);
}

#endif

