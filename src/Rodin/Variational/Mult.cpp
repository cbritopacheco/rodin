/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mult.h"

namespace Rodin::Variational
{
   Mult<ScalarFunctionBase, ScalarFunctionBase>
   operator*(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarFunctionBase, MatrixCoefficientBase>
   operator*(const ScalarFunctionBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarFunctionBase, MatrixCoefficientBase>
   operator*(const MatrixCoefficientBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(rhs, lhs);
   }

   Mult<ScalarFunctionBase, VectorCoefficientBase>
   operator*(const ScalarFunctionBase& lhs, const VectorCoefficientBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarFunctionBase, VectorCoefficientBase>
   operator*(const VectorCoefficientBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(rhs, lhs);
   }
}
