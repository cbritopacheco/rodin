#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATORUNARYMINUS_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATORUNARYMINUS_H

#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"
#include "BilinearFormIntegratorSum.h"

namespace Rodin::Variational::FormLanguage
{
   template <>
   class BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorSum>
      : public BilinearFormIntegratorSum
   {
      public:
         BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorSum& lfi);

         BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorUnaryMinus& other);

         BilinearFormIntegratorUnaryMinus(BilinearFormIntegratorUnaryMinus&&) = default;

         BilinearFormIntegratorUnaryMinus* copy() const noexcept override
         {
            return new BilinearFormIntegratorUnaryMinus(*this);
         }

      private:
         std::unique_ptr<BilinearFormIntegratorSum> m_lfi;
   };

   template <>
   class BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase>
      : public BilinearFormIntegratorBase
   {
      public:
         BilinearFormIntegratorUnaryMinus(const BilinearFormIntegratorBase& lfi);

         BilinearFormIntegratorUnaryMinus(
               const BilinearFormIntegratorUnaryMinus& other);

         const BilinearFormIntegratorBase& getBFI() const
         {
            assert(m_bfi);
            return *m_bfi;
         }

         IntegratorRegion getIntegratorRegion() const override
         {
            return getBFI().getIntegratorRegion();
         }

         const std::set<int>& getAttributes() const override
         {
            return getBFI().getAttributes();
         }

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            m_bfi->getElementMatrix(trial, test, trans, mat);
            mat *= -1.0;
         }

         BilinearFormIntegratorUnaryMinus* copy() const noexcept override
         {
            return new BilinearFormIntegratorUnaryMinus(*this);
         }

      private:
         std::unique_ptr<BilinearFormIntegratorBase> m_bfi;
   };

   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorBase> operator-(
         const BilinearFormIntegratorBase& bfi);

   BilinearFormIntegratorUnaryMinus<BilinearFormIntegratorSum> operator-(
         const BilinearFormIntegratorSum& lfi);
}

#endif

