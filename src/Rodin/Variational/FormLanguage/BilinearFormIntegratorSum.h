#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATORSUM_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATORSUM_H

#include <memory>
#include <utility>

#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   namespace Internal
   {
      class BilinearFormIntegratorSum : public mfem::BilinearFormIntegrator
      {
         public:
            BilinearFormIntegratorSum(
                  mfem::BilinearFormIntegrator& a, mfem::BilinearFormIntegrator& b)
               : m_a(a), m_b(b)
            {}

            void AssembleElementMatrix(
                  const mfem::FiniteElement& el,
                  mfem::ElementTransformation& Trans,
                  mfem::DenseMatrix& elmat)
            {
               m_a.AssembleElementMatrix(el, Trans, elmat);
               mfem::DenseMatrix tmp;
               m_b.AssembleElementMatrix(el, Trans, tmp);
               elmat.Add(1.0, tmp);
            }

         private:
            mfem::BilinearFormIntegrator& m_a;
            mfem::BilinearFormIntegrator& m_b;
      };
   }

   class BilinearFormIntegratorSum : public BilinearFormIntegratorBase
   {
      public:
         BilinearFormIntegratorSum(
               const BilinearFormIntegratorBase& lhs,
               const BilinearFormIntegratorBase& rhs);

         BilinearFormIntegratorSum(const BilinearFormIntegratorSum& other);

         BilinearFormIntegratorBase& getLHS();

         const BilinearFormIntegratorBase& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         BilinearFormIntegratorBase& getRHS();

         const BilinearFormIntegratorBase& getRHS() const
         {
            assert(m_rhs);
            return *m_rhs;
         }

         const std::set<int>& getAttributes() const override
         {
            return getLHS().getAttributes();
         }

         void buildMFEMBilinearFormIntegrator() override;

         mfem::BilinearFormIntegrator& getMFEMBilinearFormIntegrator() override;

         mfem::BilinearFormIntegrator* releaseMFEMBilinearFormIntegrator() override;

         virtual BilinearFormIntegratorSum* copy() const noexcept override
         {
            return new BilinearFormIntegratorSum(*this);
         }

      private:
         std::unique_ptr<BilinearFormIntegratorBase> m_lhs;
         std::unique_ptr<BilinearFormIntegratorBase> m_rhs;
         std::unique_ptr<Internal::BilinearFormIntegratorSum> m_mfemBFI;
   };

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs);

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs);
}

#endif

