#ifndef RODIN_VARIATIONAL_DIFFUSIONINTEGRATOR_H
#define RODIN_VARIATIONAL_DIFFUSIONINTEGRATOR_H

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   class DiffusionIntegrator : public BilinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Creates an DiffusionIntegrator with the given Lamé
          * coefficients @f$ \lambda @f$ and @f$ \mu @f$.
          */
         DiffusionIntegrator()
         {}

         DiffusionIntegrator(const DiffusionIntegrator& other)
            : BilinearFormDomainIntegrator(other),
              m_bfi(other.m_bfi)
         {}

         DiffusionIntegrator(DiffusionIntegrator&& other)
            : BilinearFormDomainIntegrator(std::move(other)),
              m_bfi(std::move(other.m_bfi))
         {}

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            assert(&trial == &test);
            m_bfi.AssembleElementMatrix(trial, trans, mat);
         }

         DiffusionIntegrator* copy() const noexcept override
         {
            return new DiffusionIntegrator(*this);
         }

      private:
         mfem::DiffusionIntegrator m_bfi;
   };
}

#endif

