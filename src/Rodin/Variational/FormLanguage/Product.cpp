/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Product.h"

namespace Rodin::Variational::FormLanguage
{
   Product<TrialFunctionBase, TestFunctionBase>
   operator*(const TrialFunctionBase& lhs, const TestFunctionBase& rhs)
   {
      return Product(lhs, rhs);
   }

   Product<ScalarCoefficientBase, TestFunctionBase>
   operator*(const ScalarCoefficientBase& lhs, const TestFunctionBase& rhs)
   {
      return Product(lhs, rhs);
   }

   Product<ScalarCoefficientBase, TestFunctionBase>
   operator*(double lhs, const TestFunctionBase& rhs)
   {
      return ScalarCoefficient(lhs) * rhs;
   }

   // ---- Product<ScalarCoefficientBase, ScalarCoefficientBase> -------------
   double
   Product<ScalarCoefficientBase, ScalarCoefficientBase>
   ::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      return getLHS().getValue(trans, ip) * getRHS().getValue(trans, ip);
   }

   Product<ScalarCoefficientBase, ScalarCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Product(lhs, rhs);
   }

   // ---- Product<ScalarCoefficientBase, MatrixCoefficientBase> -------------
   int Product<ScalarCoefficientBase, MatrixCoefficientBase>::getRows() const
   {
      return getRHS().getRows();
   }

   int Product<ScalarCoefficientBase, MatrixCoefficientBase>::getColumns() const
   {
      return getRHS().getColumns();
   }

   void Product<ScalarCoefficientBase, MatrixCoefficientBase>::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      getRHS().getValue(value, trans, ip);
      value *= getLHS().getValue(trans, ip);
   }

   Product<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return Product<ScalarCoefficientBase, MatrixCoefficientBase>(lhs, rhs);
   }

   Product<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Product<ScalarCoefficientBase, MatrixCoefficientBase>(rhs, lhs);
   }
}
