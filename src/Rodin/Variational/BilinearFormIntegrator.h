#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class BilinearFormIntegrator : public mfem::BilinearFormIntegrator
      {
         public:
            BilinearFormIntegrator(const BilinearFormIntegratorBase& bfi);

            void AssembleElementMatrix(
                  const mfem::FiniteElement& fe,
                  mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override;

            void AssembleElementMatrix2(
                  const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
                  mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override;
         private:
            std::unique_ptr<BilinearFormIntegratorBase> m_bfi;
      };
   }

   class BilinearFormIntegratorBase
      : public FormLanguage::Buildable<Internal::BilinearFormIntegrator>
   {
      public:
         virtual ~BilinearFormIntegratorBase() = default;

         /**
          * @brief Gets the attributes of the elements being integrated.
          */
         virtual const std::set<int>& getAttributes() const = 0;

         /**
          * @brief Gets the integration region.
          */
         virtual IntegratorRegion getIntegratorRegion() const = 0;

         virtual void getElementMatrix(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
         {
            assert(false);
         }


         virtual void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
         {
            assert(false);
         }

         std::unique_ptr<Internal::BilinearFormIntegrator> build() const override
         {
            return std::make_unique<Internal::BilinearFormIntegrator>(*this);
         }

         virtual BilinearFormIntegratorBase* copy() const noexcept override = 0;
   };

   class BilinearFormDomainIntegrator : public BilinearFormIntegratorBase
   {
      public:
         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         /**
          * @brief Specifies the material reference over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material reference over which the integration should
          * take place.
          */
         virtual BilinearFormDomainIntegrator& over(int attr) = 0;

         /**
          * @brief Specifies the material references over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material references over which the integration should
          * take place.
          */
         virtual BilinearFormDomainIntegrator& over(const std::set<int>& attrs) = 0;

         virtual BilinearFormDomainIntegrator* copy() const noexcept override = 0;
   };
}

#endif
