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
   /**
    * @brief Represents the Jacobian matrix @f$ \mathbf{J}_u @f$ of the
    * function @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^s \rightarrow \mathbb{R}^d @f$, the Jacobian matrix
    * @f$ \mathbf{J}_u(x) @f$ at any point @f$ x = (x_1, \ldots, x_s) @f$ is
    * defined by the @f$ s \times d @f$ matrix
    * @f[
    * \mathbf{J}_u = \begin{bmatrix}
    * \dfrac{\partial u_1}{\partial x_1} & \ldots & \dfrac{\partial u_d}{\partial x_1}\\
    * \vdots & \ddots & \vdots\\
    * \dfrac{\partial u_1}{\partial x_s} & \ldots & \dfrac{\partial u_d}{\partial x_s}
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
         constexpr
         Jacobian(GridFunction<H1>& u)
            :  m_u(u)
         {}

         constexpr
         Jacobian(const Jacobian& other)
            : m_u(other.m_u)
         {}

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint&) override
         {
            m_u.getHandle().GetVectorGradient(trans, value);
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }

      private:
         GridFunction<H1>& m_u;
   };
}

#endif
