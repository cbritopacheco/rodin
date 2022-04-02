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

   Mult<ScalarFunctionBase, MatrixFunctionBase>
   operator*(const ScalarFunctionBase& lhs, const MatrixFunctionBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarFunctionBase, MatrixFunctionBase>
   operator*(const MatrixFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(rhs, lhs);
   }

   Mult<ScalarFunctionBase, VectorFunctionBase>
   operator*(const ScalarFunctionBase& lhs, const VectorFunctionBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarFunctionBase, VectorFunctionBase>
   operator*(const VectorFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(rhs, lhs);
   }
}
