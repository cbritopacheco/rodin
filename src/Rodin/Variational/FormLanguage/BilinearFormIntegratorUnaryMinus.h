#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATORUNARYMINUS_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATORUNARYMINUS_H

#include "Rodin/Variational/BilinearFormIntegrator.h"

namespace Rodin::Variational::FormLanguage
{
   namespace Internal
   {
      class BFIUnaryMinus : public mfem::BilinearFormIntegrator
      {
         public:
            BFIUnaryMinus(mfem::BilinearFormIntegrator& bfi)
               : m_bfi(bfi)
            {}

            void AssembleElementMatrix(
                  const mfem::FiniteElement& el,
                  mfem::ElementTransformation& Trans,
                  mfem::DenseMatrix& elmat)
            {
               m_bfi.AssembleElementMatrix(el, Trans, elmat);
               elmat.Neg();
            }

         private:
            mfem::BilinearFormIntegrator& m_bfi;
      };
   }

   class BilinearFormIntegratorUnaryMinus : public BilinearFormIntegratorBase
   {
      public:
         BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorBase& lfi);

         BilinearFormIntegratorUnaryMinus(
               const BilinearFormIntegratorUnaryMinus& other);

         BilinearFormIntegratorBase& getBFI();

         const BilinearFormIntegratorBase& getBFI() const
         {
            assert(m_bfi);
            return *m_bfi;
         }

         const std::set<int>& getAttributes() const override
         {
            return getBFI().getAttributes();
         }

         void buildMFEMBilinearFormIntegrator() override;

         mfem::BilinearFormIntegrator& getMFEMBilinearFormIntegrator() override;

         mfem::BilinearFormIntegrator* releaseMFEMBilinearFormIntegrator() override;

         BilinearFormIntegratorUnaryMinus* copy() const noexcept override
         {
            return new BilinearFormIntegratorUnaryMinus(*this);
         }

      private:
         std::unique_ptr<BilinearFormIntegratorBase> m_bfi;
         std::unique_ptr<Internal::BFIUnaryMinus> m_mfemBFI;
   };

   BilinearFormIntegratorUnaryMinus operator-(
         const BilinearFormIntegratorBase& bfi);
}

#endif

