/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      /**
       * @internal
       */
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

   /**
    * @brief Represents the Jacobian matrix @f$ \mathbf{J}_u @f$ of the
    * function @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^n \rightarrow \mathbb{R}^m @f$, the Jacobian matrix
    * @f$ \mathbf{J}_u(x) @f$ at any point @f$ x = (x_1, \ldots, x_n) @f$ is
    * defined by the @f$ m \times n @f$ matrix
    * @f[
    * \mathbf{J}_u = \begin{bmatrix}
    * \dfrac{\partial u_1}{\partial x_1} & \ldots & \dfrac{\partial u_1}{\partial x_n}\\
    * \vdots & \ddots & \vdots\\
    * \dfrac{\partial u_m}{\partial x_1} & \ldots & \dfrac{\partial u_m}{\partial x_n}
    * \end{bmatrix}
    * @f]
    * 
    */
   class Jacobian : public MatrixCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Jacobian matrix of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Jacobian(GridFunction<H1>& u);

         Jacobian(const Jacobian& other);

         int getRows() const override;
         int getColumns() const override;
         void buildMFEMMatrixCoefficient() override;
         mfem::MatrixCoefficient& getMFEMMatrixCoefficient() override;

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }

      private:
         GridFunction<H1>& m_u;
         std::optional<Internal::VectorGradientCoefficient> m_mfemMatrixCoefficient;
   };
}

#endif