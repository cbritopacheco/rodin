/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Sum.h"

namespace Rodin::Variational::FormLanguage
{
   double
   Sum<ScalarCoefficientBase, ScalarCoefficientBase>
   ::getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      return getLHS().getValue(trans, ip) + getRHS().getValue(trans, ip);
   }

   void
   Sum<MatrixCoefficientBase, MatrixCoefficientBase>
   ::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      mfem::DenseMatrix m;
      m_lhs->getValue(m, trans, ip);
      m_rhs->getValue(value, trans, ip);
      value += m;
   }

   int Sum<MatrixCoefficientBase, MatrixCoefficientBase>::getRows() const
   {
      return m_lhs->getRows();
   }

   int Sum<MatrixCoefficientBase, MatrixCoefficientBase>::getColumns() const
   {
      return m_lhs->getColumns();
   }

   Sum<ScalarCoefficientBase, ScalarCoefficientBase>
   operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Sum(lhs, rhs);
   }

   Sum<MatrixCoefficientBase, MatrixCoefficientBase>
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return Sum(lhs, rhs);
   }
}
