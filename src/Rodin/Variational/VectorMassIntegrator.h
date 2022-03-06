/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORMASSINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORMASSINTEGRATOR_H

#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Omega} \lambda u \cdot v \ dx
    * @f]
    * where @f$ \lambda @f$ is a scalar coefficient.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ \lambda @f$      | ScalarCoefficient                            |
    *
    */
   class VectorMassIntegrator : public BilinearFormDomainIntegrator
   {
      public:
         /**
          * @brief Constructs a VectorMassIntegrator with a scalar coefficient
          * @f$ \lambda = 1 @f$.
          */
         VectorMassIntegrator();

         /**
          * @brief Constructs a VectorMassIntegrator with scalar coefficient
          * @f$ \lambda @f$.
          *
          * @param[in] lambda Coefficient to integrate.
          */
         VectorMassIntegrator(const ScalarCoefficientBase& lambda);

         VectorMassIntegrator(const VectorMassIntegrator& other);

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            m_mfemBFI.AssembleElementMatrix2(trial, test, trans, mat);
         }

         VectorMassIntegrator* copy() const noexcept override
         {
            return new VectorMassIntegrator(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lambda;

         std::unique_ptr<mfem::Coefficient> m_mfemLambda;
         mfem::VectorMassIntegrator m_mfemBFI;
   };

}

#endif


