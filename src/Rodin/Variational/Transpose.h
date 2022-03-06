/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRANSPOSE_H
#define RODIN_VARIATIONAL_TRANSPOSE_H

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the transpose matrix @f$ A^T @f$ of some matrix @f$ A @f$.
    *
    * For some @f$ n \times m @f$ matrix @f$ A @f$, the transpose matrix @f$
    * A^T @f$ is an @f$ m \times n @f$ matrix defined by
    * @f[
    *    {A^T}_{ij} = A_{ji}
    * @f]
    */
   class Transpose : public MatrixCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Transpose matrix of the given matrix.
          */
         Transpose(const MatrixCoefficientBase& m)
            : m_matrix(m.copy())
         {}

         Transpose(const Transpose& other)
            :  m_matrix(other.m_matrix->copy())
         {}

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            m_matrix->getValue(value, trans, ip);
            value.Transpose();
         }

         Transpose* copy() const noexcept override
         {
            return new Transpose(*this);
         }

      private:
         std::unique_ptr<MatrixCoefficientBase> m_matrix;
   };
}

#endif
