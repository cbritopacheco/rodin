/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Product.h"

namespace Rodin::Variational::FormLanguage
{
   Product<ScalarCoefficientBase, ScalarCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Product(lhs, rhs);
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

   Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>
   operator*(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
   {
      return Product(lhs, rhs);
   }
}
