/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_MATRIXCOEFFICIENTSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_MATRIXCOEFFICIENTSUM_H

#include <memory>
#include <utility>
#include <optional>

#include "Rodin/Variational/MatrixCoefficient.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class MatrixSum : public MatrixCoefficientBase
   {
      public:
         MatrixSum(
               const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs);

         MatrixSum(const MatrixSum& other);

         MatrixCoefficientBase& getLHS();
         MatrixCoefficientBase& getRHS();

         int getRows() const override;

         int getColumns() const override;

         void buildMFEMMatrixCoefficient() override;

         mfem::MatrixCoefficient& getMFEMMatrixCoefficient() override;

         virtual MatrixSum* copy() const noexcept override
         {
            return new MatrixSum(*this);
         }

      private:
         std::unique_ptr<MatrixCoefficientBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;

         std::optional<mfem::MatrixSumCoefficient> m_mfemMatrixCoefficient;
   };

   MatrixSum
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs);
}

#endif

