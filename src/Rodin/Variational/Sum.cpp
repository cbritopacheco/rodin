/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Sum.h"

namespace Rodin::Variational
{
   double
   Sum<ScalarFunctionBase, ScalarFunctionBase>
   ::getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      return getLHS().getValue(trans, ip) + getRHS().getValue(trans, ip);
   }

   void
   Sum<MatrixFunctionBase, MatrixFunctionBase>
   ::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const
   {
      mfem::DenseMatrix m;
      m_lhs->getValue(m, trans, ip);
      m_rhs->getValue(value, trans, ip);
   }

   int Sum<MatrixFunctionBase, MatrixFunctionBase>::getRows() const
   {
      return m_lhs->getRows();
   }

   int Sum<MatrixFunctionBase, MatrixFunctionBase>::getColumns() const
   {
      return m_lhs->getColumns();
   }

   Sum<ScalarFunctionBase, ScalarFunctionBase>
   operator+(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Sum(lhs, rhs);
   }

   Sum<VectorFunctionBase, VectorFunctionBase>
   operator+(const VectorFunctionBase& lhs, const VectorFunctionBase& rhs)
   {
      return Sum(lhs, rhs);
   }

   Sum<MatrixFunctionBase, MatrixFunctionBase>
   operator+(const MatrixFunctionBase& lhs, const MatrixFunctionBase& rhs)
   {
      return Sum(lhs, rhs);
   }
}
