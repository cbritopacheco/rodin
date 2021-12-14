/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class VectorGradientCoefficient : public mfem::MatrixCoefficient
      {
         public:
            VectorGradientCoefficient(mfem::GridFunction& u)
               :  mfem::MatrixCoefficient(
                     u.FESpace()->GetVDim(), u.FESpace()->GetMesh()->Dimension()),
                  m_u(u)
            {}

            virtual void Eval(
                  mfem::DenseMatrix& grad,
                  mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
            {
               T.SetIntPoint(&ip);
               m_u.GetVectorGradient(T, grad);
            }

         private:
            mfem::GridFunction& m_u;
      };
   }

   class Gradient : public MatrixCoefficientBase
   {
      public:
         Gradient(GridFunction<H1>& u);

         Gradient(const Gradient& other);

         int getRows() const override;
         int getColumns() const override;
         void buildMFEMMatrixCoefficient() override;
         mfem::MatrixCoefficient& getMFEMMatrixCoefficient() override;

         Gradient* copy() const noexcept override
         {
            return new Gradient(*this);
         }

      private:
         GridFunction<H1>& m_u;
         std::optional<Internal::VectorGradientCoefficient> m_mfemMatrixCoefficient;
   };
}

#endif
