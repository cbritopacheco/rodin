/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATORUNARYMINUS_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATORUNARYMINUS_H

#include "Rodin/Variational/LinearFormIntegrator.h"

#include "ForwardDecls.h"
#include "LinearFormIntegratorSum.h"

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
   template <>
   class LinearFormIntegratorUnaryMinus<LinearFormIntegratorSum>
      : public LinearFormIntegratorSum
   {
      public:
         LinearFormIntegratorUnaryMinus(const LinearFormIntegratorSum& lfi);

         LinearFormIntegratorUnaryMinus(const LinearFormIntegratorUnaryMinus& other);

         LinearFormIntegratorUnaryMinus(LinearFormIntegratorUnaryMinus&&) = default;

         LinearFormIntegratorUnaryMinus* copy() const noexcept override
         {
            return new LinearFormIntegratorUnaryMinus(*this);
         }

      private:
         std::unique_ptr<LinearFormIntegratorSum> m_lfi;
   };

   /**
    * @internal
    */
   template <>
   class LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
      : public LinearFormIntegratorBase
   {
      public:
         LinearFormIntegratorUnaryMinus(const LinearFormIntegratorBase& lfi);

         LinearFormIntegratorUnaryMinus(const LinearFormIntegratorUnaryMinus& other);

         LinearFormIntegratorBase& getLFI();

         const LinearFormIntegratorBase& getLFI() const
         {
            assert(m_lfi);
            return *m_lfi;
         }

         const std::set<int>& getAttributes() const override
         {
            return getLFI().getAttributes();
         }

         IntegratorRegion getIntegratorRegion() const override
         {
            return getLFI().getIntegratorRegion();
         }

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

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorBase>
   operator-(const LinearFormIntegratorBase& lfi);

   LinearFormIntegratorUnaryMinus<LinearFormIntegratorSum>
   operator-(const LinearFormIntegratorSum& lfi);
}

#endif

