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
   class Transpose : public MatrixCoefficientBase
   {
      public:
         Transpose(const MatrixCoefficientBase& m);

         Transpose(const Transpose& other);

         int getRows() const override;

         int getColumns() const override;

         void buildMFEMMatrixCoefficient() override;

         mfem::MatrixCoefficient& getMFEMMatrixCoefficient() override;

         Transpose* copy() const noexcept override
         {
            return new Transpose(*this);
         }

      private:
         std::unique_ptr<MatrixCoefficientBase> m_matrix;
         std::optional<mfem::TransposeMatrixCoefficient> m_mfemMatrixCoefficient;
   };
}

#endif
